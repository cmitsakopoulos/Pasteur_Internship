import os
import sys
import datetime
import pandas as pd
from abc import ABC, abstractmethod
from typing import Dict

from function_dump import (
    extractor,
    parser,
    calculate_cdr_chars,
    calculate_antigen_chars,
)

#Separation of all application functions into steps, which are brought together in a recipe like way, to accomplish different tasks; aka chem, clean, or extract from version 1 (=version lame) of this program.
#Object oriented approach, where dfs are the key object being traded between classes/functions, enables more accesible operation with different file formats. I cannot train my database solely with a former parser that reads only NAstructural database csvs.

class Pipeline:
    def __init__(self, steps: list[Step]):
        self.steps = steps
        self.data: Dict[str, pd.DataFrame] = {}

    def run(self) -> Dict[str, pd.DataFrame]:
        for step in self.steps:
            out = step.process(self.data)
            # merge new entries into the central data dict
            self.data.update(out)
        return self.data
    
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

        return result
    
class Parser(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        # Identify keys corresponding to raw CSV inputs (ie csv1, csv2, etc.)
        csv_keys = [key for key in data if key.startswith('csv')]
        for key in csv_keys:
            pre_parsed_df = data.pop(key)
            antigen_df, cdr_df = parser(pre_parsed_df)
            new_data[f"antigen_{key}"] = antigen_df
            new_data[f"cdr_{key}"] = cdr_df
        return new_data


class Write(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        # Ensure directory exists; else we die
        output_dir = os.path.join(os.path.dirname(__file__), "Computation_Deposit")
        os.makedirs(output_dir, exist_ok=True)
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        for key, df in data.items():
            filename = f"{key}_{ts}.csv"
            filepath = os.path.join(output_dir, filename)
            df.to_csv(filepath, index=False)
            print(f"Wrote DataFrame '{key}' to {filepath}")
        return {}

#PLACEHOLDER AI GENERATED, NEED TO THINK ABOUT THIS FUNCTIONALITY
class InspectionStep(Step):
    pass

#CHEMICAL CHARACTERISTICS
class CDRComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'cdr' in data:
            df = data['cdr']
            calculate_cdr_chars(df)
            new_data['cdr'] = df
        return new_data

class AntigenComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'antigen' in data:
            df = data['antigen']
            calculate_antigen_chars(df)
            new_data['antigen'] = df
        return new_data
#CHEMICAL CHARACTERISTICS

# Multiple antibodies can bind to a single antigen, many antigens can sature a single antibody, this is reflected in the source data too...
class FlattenDuplicates(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'antigen' in data:
            df = data['antigen']
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
        return new_data

#Named after the os import function "walk", traverses user provided directory to build our shared dict, which is then parsed further...
class Walker(Step):
    def __init__(self, input_path: str):
        self.input_path = input_path
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        count = 1
        for root, _, files in os.walk(self.input_path):
            for fname in files:
                if not fname.lower().endswith('.csv'):
                    continue
                filepath = os.path.join(root, fname)
                try:
                    df = extractor(filepath)
                except Exception as e:
                    print(f"Warning: failed to extract {filepath!r}: {e}", file=sys.stderr)
                    continue
                key = f"csv{count}"
                new_data[key] = df
                count += 1
        return new_data

#FUNCTIONS; keep in mind that due to the idiosyncratic formatting of NAStructuralDB, we need to increase max memory usage for csvs....I made it unlimited, change if needed
