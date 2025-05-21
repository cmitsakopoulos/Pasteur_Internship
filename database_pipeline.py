import os
import datetime
import pandas as pd
import re
from abc import ABC, abstractmethod
from typing import Dict
from concurrent.futures import ProcessPoolExecutor

from function_dump import (
    extractor,
    parser,
    calculate_cdr_chars,
    calculate_antigen_chars,
)

#Separation of all application functions into steps, which are brought together in a recipe like way, to accomplish different tasks; aka chem, clean, or extract from version 1 (=version lame) of this program.
#Object oriented approach, where dfs are the key object being traded between classes/functions, enables more accesible operation with different file formats. I cannot train my database solely with a former parser that reads only NAstructural database csvs.
#CLI is connected and feeds a list of steps depending on the type of recipe the user wants to use: ex. complete to do a full treatment of teh raw data, or simple for testing purposes. Inspect clauses are handled independently. I preferred inspection with manually created functions, because pandas is nicer than dunder methods. 

class Pipeline: #Receives a "recipe"/list of Step subclasses needed to produce a desired output
    def __init__(self, steps):
        self.steps = steps

    def run(self) -> Dict[str, pd.DataFrame]:
        data: Dict[str, pd.DataFrame] = {}
        for step in self.steps:
            data = step.process(data)
        return data
    
class Step(ABC): #Parent class to be method overrriden by specialised children
    @abstractmethod
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        pass

class Concatenation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        cdr_frames = []
        antigen_frames = []
        for key, df in data.items():
            if key.startswith('cdr_'):
                cdr_frames.append(df)
            elif key.startswith('antigen_'):
                antigen_frames.append(df)
        result: Dict[str, pd.DataFrame] = {}
        if cdr_frames:
            result['cdr'] = pd.concat(cdr_frames, ignore_index=True)
        else:
            result['cdr'] = pd.DataFrame()
            print("No CDRS, moving on")
        if antigen_frames:
            result['antigen'] = pd.concat(antigen_frames, ignore_index=True)
        else:
            result['antigen'] = pd.DataFrame()
            print("No Antigens, moving on")

        print(f"Concatenation → output keys: {list(result.keys())}")
        return result

class Write(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        # Ensure directory exists; else we die
        output_dir = os.path.join(os.path.dirname(__file__), "Computation_Deposit")
        os.makedirs(output_dir, exist_ok=True)
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") #can remove if you want, obviously
        for key, df in data.items():
            filename = f"{key}_{ts}.csv"
            filepath = os.path.join(output_dir, filename)
            df.to_csv(filepath, index=False)
            print(f"Wrote DataFrame '{key}' to {filepath}")
        print(f"Write → wrote {len(data)} tables: {list(data.keys())}")
        return {}

#CHEMICAL CHARACTERISTICS
class CDRComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = data.copy() #Copy entire thing then mutate
        #Catch the antigen dataframes, remember that Walker will name the files antigen_csv or cdr_csv something
        if 'cdr' in new_data:
            computation_cdr_df = new_data["cdr"]
            calculate_cdr_chars(computation_cdr_df)
            new_data['cdr'] = computation_cdr_df
        if 'cdr' in new_data:
            print(f"CDRComputation → processed cdr, rows: {new_data['cdr'].shape[0]}")
        return new_data #return "mutated" dataframe, feed back into the system

class AntigenComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = data.copy() 
        if 'antigen' in data:
            computation_antigen_df = new_data["antigen"]
            calculate_antigen_chars(computation_antigen_df)
            new_data['antigen'] = computation_antigen_df
        if 'antigen' in new_data:
            print(f"AntigenComputation → processed antigen, rows: {new_data['antigen'].shape[0]}")
        return new_data
#CHEMICAL CHARACTERISTICS

# Multiple antibodies can bind to a single antigen, many antigens can sature a single antibody, this is reflected in the source data too...
class FlattenDuplicates(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'cdr' in data:
            cdr_df = data['cdr'].copy()
            if not cdr_df['h3_chain'].duplicated().any(): #Forgot to add a clause by which if there are no duplicates, we skip...AND return the df as is
                print(f"FlattenDuplicates → no CDR duplicates found, keeping {len(cdr_df)} rows")
                new_data['cdr'] = cdr_df
            else:
                dict = {
                    col: 'first'
                    for col in cdr_df.columns
                    if col not in ('h3_chain', 'pdb_id')
                }
                dict['pdb_id'] = lambda ids: list(ids.unique())
                flat_cdr = (
                    cdr_df
                    .groupby('h3_chain')
                    .agg(dict)
                    .reset_index()
                )
                new_data['cdr'] = flat_cdr
        if 'antigen' in data:
            df = data['antigen'].copy()
            if not df['antigen_seq'].duplicated().any(): 
                print(f"FlattenDuplicates → no Antigen duplicates found, keeping {len(df)} rows")
                new_data['antigen'] = df
            else:
                agg_dict = {
                    col: 'first'
                    for col in df.columns
                    if col not in ('antigen_seq', 'antigen_id')
                }
                agg_dict['antigen_id'] = lambda ids: list(ids.unique())
                flat_df = (
                    df
                    .groupby('antigen_seq')
                    .agg(agg_dict)
                    .reset_index()
                )
                new_data['antigen'] = flat_df
        if 'antigen' in new_data:
            print(f"FlattenDuplicates → flattened antigen, rows: {new_data['antigen'].shape[0]}")
        if 'cdr' in new_data:
            print(f"FlattenDuplicates → flattened CDR, rows: {new_data['antigen'].shape[0]}")
        return new_data
    
#Named after the os import function "walk", traverses user provided directory to build our shared dict, which is then parsed further...
#Parallelisation method to improve computational time is AI assisted code, any problems that arise are not my fault; none discovered for now.
class Walker(Step):
    def __init__(self, input_path: str):
        self.input_path = input_path

    @staticmethod
    def _process_file(args):
        filepath, idx = args
        df_raw = extractor(filepath)
        antigen_df, cdr_df = parser(df_raw)
        return filepath, idx, antigen_df, cdr_df

    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        tasks = []
        count = 0
        for root, _, files in os.walk(self.input_path):
            for fname in files:
                if not fname.lower().endswith('.csv'):
                    continue
                count += 1
                filepath = os.path.join(root, fname)
                tasks.append((filepath, count))
        with ProcessPoolExecutor(max_workers=2) as executor:
            for filepath, idx, antigen_df, cdr_df in executor.map(self._process_file, tasks):
                base_key = f"csv{idx}"
                new_data[f"antigen_{base_key}"] = antigen_df
                new_data[f"cdr_{base_key}"]      = cdr_df
                print(
                    f"Walker → processed {filepath!r}: "
                    f"antigen rows={antigen_df.shape[0]}, "
                    f"cdr rows={cdr_df.shape[0]}"
                )
        return new_data

#Offering a simple to implement solution to a complicated issue: remove purification tags from antigenic sequences no matter if they are alone or in groups or whatever.
class RmPurificationTags(Step):
    def __init__(self):#Each class is only instantiated once as it is a componenent in the larger recipe called by the cli, so this, which is normally bad practice, ill leave, as it makes no difference due to being instantiated just once/just how I know how to do it...
        # we compile the pattern for motifs once, then it looks for it hopefully finds it..
        self.motifs = {
            "c-Myc": "EQKLISEEDL",
            "6-His": "HHHHHH",
            "5-His": "HHHHH", #Seems like in the NAStructural DB, some of the antigenic sequences have N or C-terminal pentameric His-tags attached...accidental clipping of the sequence at the ends?!
            "FLAG":  "DYKDDDDK",
            "V5":    "GKPIPNPLLGLDST",
        }
        self.max_gap = 6 #Change if you want to make the gap between purifiaction tag motifs smaller...
        #Preparing for alteration (aka so that "re" can look for c-Myc OR '|' FLAG OR/AND 6-His)
        self.alteration = "|".join(map(re.escape, self.motifs.values()))
        # one motif, followed by 0-N chunks of ≤7 AAs + another motif
        self.clusters = re.compile(
            rf'(?:{self.alteration}(?:[A-Z]{{1,{self.max_gap}}}{self.alteration})*)'
        ) #create a rubric for the purging to occur later
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if "cdr" in data: # pass cdr through unchanged
            new_data["cdr"] = data["cdr"].copy()
        if "antigen" in data:
            df = data["antigen"].copy()
            for idx, seq in df["antigen_seq"].items():
                df.at[idx, "antigen_seq"] = self.clusters.sub("", seq)
            new_data["antigen"] = df
        return new_data
    
class AssignIDs(Step):
    def __init__(self): #Each class is only instantiated once as it is a componenent in the larger recipe called by the cli, so this which is normally bad practice, ill leave as it makes no difference due to being instantiated just once/just how I know how to do it...
        self.amino_acid_rubric = {
        'A': 1,  'C': 2,  'D': 3,  'E': 4,  'F': 5,
        'G': 6,  'H': 7,  'I': 8,  'K': 9,  'L':10,
        'M':11,  'N':12,  'P':13,  'Q':14,  'R':15,
        'S':16,  'T':17,  'V':18,  'W':19,  'Y':20,}
        self.re_pattern = re.compile(r'[ACDEFGHIKLMNPQRSTVWY]')
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'cdr' in data:
            cdr_df = data["cdr"].copy()
            cdr_df["cdr_computed_id"] = pd.NA #Introduce new column for our computed id
            for idx, seq in cdr_df["h3_chain"].items():
                if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq): #Check if the filtering done before didnt work; maybe we dont only have ambiguous X, maybe there are additional characters
                    raise ValueError(f"Illegal residue in {seq!r} at row {idx}")
                else:
                    numeric_str = self.re_pattern.sub(lambda m: f"{self.amino_acid_rubric[m.group(0)]}+", seq) #Replace AAs with encodings
                    cdr_df.at[idx, "cdr_computed_id"] = sum(map(int, numeric_str.rstrip("+").split("+"))) #Sum the AA encodings
            new_data["cdr"] = cdr_df
        if "antigen" in data:
            ant_df = data["antigen"].copy()
            ant_df["antigen_computed_id"] = pd.NA 
            for idx, seq in ant_df["antigen_seq"].items():
                if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq):
                    raise ValueError(f"Illegal residue in {seq!r} at row {idx}")
                else:
                    numeric_str = self.re_pattern.sub(lambda m: f"{self.amino_acid_rubric[m.group(0)]}+", seq)
                    ant_df.at[idx, "antigen_computed_id"] = sum(map(int, numeric_str.rstrip("+").split("+")))
            new_data["antigen"] = ant_df
        return new_data