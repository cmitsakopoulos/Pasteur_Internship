import os
import datetime
from concurrent.futures import ProcessPoolExecutor
import re
from abc import ABC, abstractmethod
from typing import Dict
import pandas as pd
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy import text
from sqlalchemy.pool import NullPool



from function_dump import (
    extractor,
    parser,
    calculate_cdr_chars,
    calculate_antigen_chars,
    to_list,
    prefix,
    clear_dir,
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
        output_dir = os.path.join(os.path.dirname(__file__), "Computation_Deposit")
        os.makedirs(output_dir, exist_ok=True)
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
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
        new_data: Dict[str, pd.DataFrame] = data.copy()
        if 'cdr' in new_data:
            computation_cdr_df = new_data["cdr"]
            calculate_cdr_chars(computation_cdr_df)
            new_data['cdr'] = computation_cdr_df
        if 'cdr' in new_data:
            print(f"CDRComputation → processed cdr, rows: {new_data['cdr'].shape[0]}")
        return new_data

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
            if not cdr_df['h3_chain'].duplicated().any():
                print(f"FlattenDuplicates → no CDR duplicates found, keeping {len(cdr_df)} rows")
                new_data['cdr'] = cdr_df
            else:
                search_dict = {
                    col: 'first'
                    for col in cdr_df.columns
                    if col not in ('h3_chain', 'pdb_id')
                }
                search_dict['pdb_id'] = lambda ids: list(ids.unique())
                flat_cdr = (
                    cdr_df
                    .groupby('h3_chain')
                    .agg(search_dict)
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
                    if col not in ('antigen_seq', 'corresponding_pdb_antibody')
                }
                agg_dict['corresponding_pdb_antibody'] = lambda ids: list(ids.unique())
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
            print(f"FlattenDuplicates → flattened CDR, rows: {new_data['cdr'].shape[0]}")
        return new_data

class CleanUp(Step): #Connected to function in function dump, just cleans up the code here
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        internal_dir = os.path.join(os.path.dirname(__file__), "Internal_Files")
        output_dir = os.path.join(os.path.dirname(__file__), "Computation_Deposit")
        clear_dir(internal_dir)
        clear_dir(output_dir)
        data = {}
        return data

class PreWalker(Step): #Not incredibly efficient but very flexible method to check if files have already been parsed/extracted...it will also read files within the common directory if needed.
    def __init__(self, input_path: str):
        self.input_path = Path(input_path)
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        internal_dir = Path(__file__).parent / "Internal_Files"
        input_files = list(self.input_path.rglob("*.csv"))
        if internal_dir.is_dir():
            processed_stems = {p.stem for p in internal_dir.rglob("*.csv")}
            tasks_to_be = [
                str(f.resolve())
                for f in input_files
                if f.stem not in processed_stems
            ]
            dropped = len(list(internal_dir.glob('*.csv')))
            data["paths"] = tasks_to_be
            print(f"Dropped {dropped} already-processed files; {len(tasks_to_be)} remain.")
            for p in internal_dir.rglob("*.csv"):
                stem = p.stem.lower()
                df = pd.read_csv(p)
                if "antigen" in stem:
                    data[f"antigen_{stem}"] = df
                elif "cdr" in stem:
                    data[f"cdr_{stem}"] = df
        else:
            data["paths"] = [str(f.resolve()) for f in input_files]
        return data
    
#Named after the os import function "walk", traverses user provided directory to build our shared dict, which is then parsed further...
class Walker(Step):
    @staticmethod
    def _process_file(args):
        filepath, idx = args
        df_raw = extractor(filepath)
        antigen_df, cdr_df = parser(df_raw)
        return filepath, idx, antigen_df, cdr_df
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        tasks = []
        output_dir_pf = os.path.join(os.path.dirname(__file__), "Internal_Files")
        os.makedirs(output_dir_pf, exist_ok=True)
        for idx, filepath in enumerate(data["paths"], start=1):
            tasks.append((filepath, idx))
        with ProcessPoolExecutor(max_workers=2) as executor:
            for filepath, idx, antigen_df, cdr_df in executor.map(self._process_file, tasks):
                base_key = Path(filepath).stem
                if not antigen_df.empty:
                    new_data[f"antigen_{base_key}"] = antigen_df
                    filename_agn = f"{base_key}_antigen.csv"
                    filepath_agn = os.path.join(output_dir_pf, filename_agn)
                    antigen_df.to_csv(filepath_agn, index=False)
                    print(f"Wrote Antigen DataFrame from '{base_key}' to {filepath_agn}")
                if not cdr_df.empty:
                    new_data[f"cdr_{base_key}"] = cdr_df
                    filename_abd = f"{base_key}_cdr.csv"
                    filepath_abd = os.path.join(output_dir_pf, filename_abd)
                    cdr_df.to_csv(filepath_abd, index=False)
                    print(f"Wrote CDR DataFrame from '{base_key}' to {filepath_abd}")
                print(
                    f"Walker → processed {filepath!r}: "
                    f"antigen rows={antigen_df.shape[0]}, "
                    f"cdr rows={cdr_df.shape[0]}"
                )
        return new_data

#Offering a simple to implement solution to a complicated issue: remove purification tags from antigenic sequences no matter if they are alone or in groups or whatever.
class RmPurificationTags(Step):
    def __init__(self):
        self.motifs = {
            "c-Myc": "EQKLISEEDL",
            "6-His": "HHHHHH",
            "5-His": "HHHHH",
            "FLAG":  "DYKDDDDK",
            "V5":    "GKPIPNPLLGLDST",
        }
        self.max_gap = 6
        self.alteration = "|".join(map(re.escape, self.motifs.values()))
        self.clusters = re.compile(
            rf'(?:{self.alteration}(?:[A-Z]{{1,{self.max_gap}}}{self.alteration})*)'
        )
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if "cdr" in data:
            new_data["cdr"] = data["cdr"].copy()
        if "antigen" in data:
            df = data["antigen"].copy()
            for idx, seq in df["antigen_seq"].items():
                df.at[idx, "antigen_seq"] = self.clusters.sub("", seq)
            new_data["antigen"] = df
            print(f"RmPurificationTags → cleaned antigenic sequences in the following number of rows: {new_data['antigen'].shape[0]}")
        return new_data

#Remember, this step is ONLY intended for assigning unique identifiers, this is uninformative encoding, all prior information is lost dramatically; the sequence cannot possibly be discerned from the computed ID. 
class AssignIDs(Step):
    def __init__(self):
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
            cdr_df["cdr_computed_id"] = pd.NA
            for idx, seq in cdr_df["h3_chain"].items():
                if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq):
                    raise ValueError(f"Illegal residue in {seq!r} at row {idx}")
                else:
                    cdr_df.at[idx, "cdr_computed_id"] = self.re_pattern.sub(lambda m: str(self.amino_acid_rubric[m.group(0)]), seq)
            new_data["cdr"] = cdr_df
            print(f"AssignIDs → assigned {new_data['cdr'].shape[0]} sequence based IDs")
        if "antigen" in data:
            ant_df = data["antigen"].copy()
            ant_df["antigen_computed_id"] = pd.NA
            for idx, seq in ant_df["antigen_seq"].items():
                if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq):
                    raise ValueError(f"Illegal residue in {seq!r} at row {idx}")
                else:
                    ant_df.at[idx, "antigen_computed_id"] = self.re_pattern.sub(lambda m: str(self.amino_acid_rubric[m.group(0)]), seq)
            new_data["antigen"] = ant_df
            print(f"AssignIDs → assigned {new_data['antigen'].shape[0]} sequence based IDs")
        return new_data


#Bear in mind this algorithm prefers the use of dictionaries to compact data; might not be entirely sound..
#Importantly, one could use pandas to make use of the C based hash join?! However I am more comfortable with dictionaries and feel that an inner join with pandas might miss out on the multiplicity of connections between different antigen/cdr objects... this is entirely based on vibes
#Algorithmic complexity is therefore higher than pandas but whatever...preliminary trials showed that to parse throuh the entire NAStructuralDB deduplicated dataset and find 679 relationships, it takes 5 seconds on M1 Macintosh without connection to power supply!
class ComputeRelationships(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        antigens = data["antigen"].copy()
        cdrs = data["cdr"].copy()
        antigens["pdb_list"] = antigens["corresponding_pdb_antibody"].apply(to_list)
        cdrs["pdb_list"]     = cdrs["pdb_id"].apply(to_list)
        ag_exp  = antigens.explode("pdb_list")[["antigen_computed_id", "pdb_list"]]
        cdr_exp = cdrs.explode("pdb_list")[["cdr_computed_id",     "pdb_list"]]
        pairs = (
            pd.merge(ag_exp, cdr_exp, on="pdb_list", how="inner")
              [["antigen_computed_id", "cdr_computed_id"]]
              .drop_duplicates()
              .reset_index(drop=True)
        )
        new_data["antigen"]       = antigens.drop(columns=["pdb_list"])
        new_data["cdr"]           = cdrs.drop(columns=["pdb_list"])
        new_data["relationships"] = pairs
        return new_data
    
class WorkWithDatabase(Step):
    def __init__(self):
        self.connection_string = "postgresql://chrismitsacopoulos:password@localhost/pasteurdb"
        self.ddl_dir = Path(__file__).resolve().parent / 'sql_files'
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        ddl_paths = sorted(self.ddl_dir.glob("*.sql"), key=prefix)
        if not ddl_paths:
            raise RuntimeError(f"No DDL files found in {self.ddl_dir}")
        engine = create_engine(self.connection_string, echo=True, poolclass=NullPool, pool_pre_ping=True)
        with engine.begin() as conn:
            first = ddl_paths[0]
            print(f"WorkWithDatabase → creating staging tables from {first.name}")
            for stmt in first.read_text(encoding="utf-8").split(";"):
                stmt = stmt.strip()
                if stmt and not re.fullmatch(r"(?i)(BEGIN|COMMIT|ROLLBACK)", stmt):
                    conn.execute(text(stmt))
            print("WorkWithDatabase → inserting DataFrames into staging tables...")
            for dict_key, dict_df in data.items():
                if not dict_df.empty:
                    stg_table = f"staging_{dict_key}"
                    print(f" • loading {len(dict_df)} rows into {stg_table}")
                    dict_df.to_sql(stg_table, conn, if_exists="append", index=False)
            print("WorkWithDatabase → altering staging_antigen.antigen_is_incomplete to BOOLEAN...")
            conn.execute(text(      """
                ALTER TABLE staging_antigen
                  ALTER COLUMN antigen_is_incomplete TYPE BOOLEAN
                  USING (antigen_is_incomplete = 1);
                """)
            )
            print("WorkWithDatabase → altering staging_cdr.h3_is_incomplete to BOOLEAN...")
            conn.execute(text(         """
                ALTER TABLE staging_cdr
                  ALTER COLUMN h3_is_incomplete TYPE BOOLEAN
                  USING (h3_is_incomplete = 1);
                """)
            )
            print("WorkWithDatabase → altering staging_cdr.l3_is_incomplete to BOOLEAN...")
            conn.execute(text( """
                ALTER TABLE staging_cdr
                  ALTER COLUMN l3_is_incomplete TYPE BOOLEAN
                  USING (l3_is_incomplete = 1);
                """)
            )
            for ddl_path in ddl_paths[1:]:
                print(f"WorkWithDatabase → applying post-load DDL {ddl_path.name}")
                for stmt in ddl_path.read_text(encoding="utf-8").split(";"):
                    stmt = stmt.strip()
                    if stmt and not re.fullmatch(r"(?i)(BEGIN|COMMIT|ROLLBACK)", stmt):
                        conn.execute(text(stmt))
        print("WorkWithDatabase → all DDL applied and data loaded.")
        return data