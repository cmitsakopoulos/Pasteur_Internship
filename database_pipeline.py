import os
import sys
import datetime
import pandas as pd
import click
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

class Pipeline:
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
        new_data: Dict[str, pd.DataFrame] = {}
        if 'antigen' in data:
            new_data['antigen'] = data['antigen']
        if 'cdr' in data:
            df = data['cdr']
            calculate_cdr_chars(df)
            new_data['cdr'] = df
        if 'cdr' in new_data:
            print(f"CDRComputation → processed cdr, rows: {new_data['cdr'].shape[0]}")
        return new_data

class AntigenComputation(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'cdr' in data:
            new_data['cdr'] = data['cdr']
        if 'antigen' in data:
            df = data['antigen']
            calculate_antigen_chars(df)
            new_data['antigen'] = df
        if 'antigen' in new_data:
            print(f"AntigenComputation → processed antigen, rows: {new_data['antigen'].shape[0]}")
        return new_data
#CHEMICAL CHARACTERISTICS

# Multiple antibodies can bind to a single antigen, many antigens can sature a single antibody, this is reflected in the source data too...
class FlattenDuplicates(Step):
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        if 'cdr' in data:
            new_data['cdr'] = data['cdr']
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
        if 'antigen' in new_data:
            print(f"FlattenDuplicates → flattened antigen, rows: {new_data['antigen'].shape[0]}")
        return new_data
#Named after the os import function "walk", traverses user provided directory to build our shared dict, which is then parsed further...
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
        # 2) Parallel extraction + parsing
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