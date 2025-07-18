import os
import sys
import csv
import json
import re
import shutil
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import umap
from propy.PyPro import GetProDes


from app_components.config import INTERNAL_FILES_DIR

#For visual purposes, Ive split the functions into a separate file so that it is easier to debug...only personal, literally-obviously zero difference...

csv.field_size_limit(500_000_000)

def extractor(filepath):
    if not os.path.isfile(filepath):
        print(f"Error: file not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    else:
        records = []
        with open(filepath, newline='', encoding='utf-8') as f:
            reader = csv.reader(
                f,
                delimiter=',',
                quotechar='"',
                escapechar='\\',
                doublequote=False
            )
            try:
                headers = next(reader)
            except StopIteration:
                print("Error: CSV file is empty", file=sys.stderr)
                sys.exit(1)
            expected = ["json", "pdb_id", "heavy", "light", "antigen", "cl_id"]
            if headers != expected: #Fallback to normal csv reading if its not the diffcult one
                print("NORMAL csv type, reverting to standard CSV read.", file=sys.stderr)
                non_idiosyncratic_csv = pd.read_csv(filepath)
                return  non_idiosyncratic_csv #If its not idiosyncratic
            for lineno, row in enumerate(reader, start=2):
                if len(row) != 6:
                    print(f"Line {lineno}: expected 6 columns, got {len(row)}", file=sys.stderr)
                    continue
                json_blob, pdb_id, heavy, light, antigen, cl_id = row
                try:
                    parsed = json.loads(json_blob)
                except json.JSONDecodeError as e:
                    print(f"Line {lineno}: JSON error: {e}", file=sys.stderr)
                    continue
                parsed.update({
                    "pdb_id": pdb_id,
                    "heavy": heavy,
                    "light": light,
                    "antigen": antigen,
                    "cl_id": cl_id
                })
                records.append(parsed)
    return pd.DataFrame(records)

def parser(df):
    cdr_records = []
    antigen_records = []

    if 'basic' in df.columns: # Handles the idiosyncratic format
        for idx, basic_dict in enumerate(df.get('basic', [])):
            cdr_row = {}
            antigen_row = {}
            if isinstance(basic_dict, dict):
                cdr_row['pdb_id'] = basic_dict.get('pdb_id')
                cdr_row['database_origin'] = "NAStructuralDB"
                antigen_row['database_origin'] = "NAStructuralDB"
                antigen_row['antigen_seq'] = basic_dict.get('antigen_seq')
                cdr_row['heavy_taxonomy_id'] = basic_dict.get('heavy_taxonomy_id')
                cdr_row['heavy_host_organism_name'] = basic_dict.get('heavy_host_organism_name')
                antigen_row['antigen_organism_name'] = basic_dict.get('antigen_organism_name')
                antigen_row['antigen_host_organism'] = basic_dict.get('antigen_host_organism')
                antigen_row['antigen_taxonomy_id'] = basic_dict.get('antigen_taxonomy_id')
                cdr_row['resolution'] = basic_dict.get('resolution')
                cdr_row['method'] = basic_dict.get('method')
                cdr_row['last_update'] = basic_dict.get('last_update')
                antigen_row['resolution'] = basic_dict.get('resolution')
                antigen_row['method'] = basic_dict.get('method')
                antigen_row['last_update'] = basic_dict.get('last_update')
                l3_h3 = basic_dict.get('l3_h3', '')
                parts = l3_h3.split('_', 1)
                cdr_row['l3_chain'] = parts[0] if len(parts) > 0 else None
                cdr_row['h3_chain'] = parts[1] if len(parts) > 1 else None
                antigen_row['corresponding_pdb_antibody'] = f"{basic_dict.get('pdb_id')}"
            cdr_records.append(cdr_row)
            antigen_records.append(antigen_row)
    else: 
        column_map = {}
        try:
            instructions_path = INTERNAL_FILES_DIR / "csv_instructions.json"
            with open(instructions_path, 'r') as f:
                column_map = json.load(f)
        except FileNotFoundError:
            comment = f"Warning: Instructions file not found at '{instructions_path}'. Failing."
            raise FileNotFoundError(comment)
        except json.JSONDecodeError:
            comment_2 = f"Warning: Could not decode JSON from '{instructions_path}'. Failing."
            raise ValueError(comment_2) 
        for idx, new_row in df.iterrows():
            cdr_row = {}
            antigen_row = {}
            pdb_col = column_map.get('pdb_name')
            pdb_val = new_row.get(pdb_col) if pdb_col else pd.NA
            dataset_col = column_map.get('dataset')
            dataset_val = new_row.get(dataset_col) if dataset_col else "User_Upload"
            antigen_row['corresponding_pdb_antibody'] = pdb_val
            cdr_row['pdb_id'] = pdb_val
            antigen_seq_col = column_map.get('antigen_sequence')
            cdr_h3_col = column_map.get('heavy_cdr3')
            cdr_l3_col = column_map.get('light_cdr3')
            antigen_row["antigen_seq"] = new_row.get(antigen_seq_col) if antigen_seq_col else pd.NA
            cdr_row["h3_chain"] = new_row.get(cdr_h3_col) if cdr_h3_col else pd.NA
            cdr_row['l3_chain'] = new_row.get(cdr_l3_col) if cdr_l3_col else pd.NA
            cdr_row['database_origin'] = dataset_val
            antigen_row['database_origin'] = dataset_val
            cdr_row['heavy_taxonomy_id'] = pd.NA
            cdr_row['heavy_host_organism_name'] = pd.NA
            antigen_row['antigen_organism_name'] = pd.NA
            antigen_row['antigen_host_organism'] = pd.NA
            antigen_row['antigen_taxonomy_id'] = pd.NA
            cdr_row['resolution'] = pd.NA
            cdr_row['method'] = pd.NA
            cdr_row['last_update'] = pd.NA
            antigen_row['resolution'] = pd.NA
            antigen_row['method'] = pd.NA
            antigen_row['last_update'] = pd.NA
            cdr_records.append(cdr_row)
            antigen_records.append(antigen_row)
            
    antigen_df = pd.DataFrame(antigen_records)
    cdr_df = pd.DataFrame(cdr_records)
    return antigen_df, cdr_df      

def calculate_cdr_chars(df: pd.DataFrame) -> pd.DataFrame:

    df_out = df.copy()

    if 'h3_chain' in df_out.columns:
        h3_original = df_out['h3_chain'].dropna()
        h3_clean = h3_original.str.replace('X', '')
        
        df_out['h3_is_incomplete'] = h3_original.str.contains('X').astype(int)
        df_out['h3_pi'] = h3_clean.apply(lambda s: ProteinAnalysis(s).isoelectric_point())
        df_out['h3_ctd'] = h3_clean.apply(compute_ctd_vector)
        df_out['h3_n_glycosilation_sites'] = h3_original.apply(find_N_glycosilation)
        df_out['h3_o_glycosilation_sites'] = h3_original.apply(
            lambda s: len(re.findall(r'(?=(?:[ST]{3}))(?!P)', s)) if isinstance(s, str) else 0
        )
        df_out['h3_chain'] = h3_clean

    if 'l3_chain' in df_out.columns:
        l3_original = df_out['l3_chain'].dropna()
        l3_clean = l3_original.str.replace('X', '')
        df_out['l3_is_incomplete'] = l3_original.str.contains('X').astype(int)
        df_out['l3_pi'] = l3_clean.apply(lambda s: ProteinAnalysis(s).isoelectric_point())
        df_out['l3_n_glycosilation_sites'] = l3_original.apply(find_N_glycosilation)
        df_out['l3_o_glycosilation_sites'] = l3_original.apply(
            lambda s: len(re.findall(r'(?=(?:[ST]{3}))(?!P)', s)) if isinstance(s, str) else 0
        )
        df_out['l3_chain'] = l3_clean

    return df_out


# O(n) complexity algorithm, I will not concede to the re automated functions
def find_N_glycosilation(seq):
    seq_length = len(seq)
    i = 0
    count = 0
    while i+2 < seq_length:
        if seq[i] == "N" and seq[i+1] != "P" and seq[i+2] in ("T", "S"):
            count += 1
            i += 1
        else:
            i += 1
    return count

def calculate_antigen_chars(df: pd.DataFrame) -> pd.DataFrame:
    df_out = df.copy()

    if 'antigen_seq' in df_out.columns:
        antigen_original = df_out['antigen_seq'].dropna()
        antigen_clean = antigen_original.str.replace('X', '')
        df_out['antigen_is_incomplete'] = antigen_original.str.contains('X').astype(int)
        df_out['antigen_pi'] = antigen_clean.apply(lambda s: ProteinAnalysis(s).isoelectric_point())
        df_out['ant_ctd'] = antigen_clean.apply(compute_ctd_vector)
        df_out['antigen_seq'] = antigen_clean
        
    return df_out

def inspect_verbose(df):
    for idx, row in df.iterrows():
        print(f"\nRow {idx}:")
        for key, val in row.items():
            print(f"  {key}: {val}")

def inspect_summary(df):
    print(f"Provided DataFrame contains {len(df)} records")

def to_list(x): 
    if isinstance(x, (list, tuple, set, pd.Series)):
        return list(x)           
    if pd.isna(x):
        return []
    return [x]

def prefix(path): 
    m = re.match(r"(\d+)", path.stem)
    if not m:
        raise ValueError(f"SQL filename '{path.name}' missing numeric prefix")
    return int(m.group(1))

def clear_dir(path: str, exception=None):
    if not os.path.isdir(path):
        return
    for name in os.listdir(path):
        if name == exception:
            continue
        full = os.path.join(path, name)
        if os.path.isfile(full) or os.path.islink(full):
            os.remove(full)
        elif os.path.isdir(full):
            shutil.rmtree(full)

def UMAP(distance_matrix: np.ndarray) -> np.ndarray:
    reducer = umap.UMAP(
        n_components=2,
        metric='precomputed',
        random_state=42
    )
    coords_2d = reducer.fit_transform(distance_matrix)
    return coords_2d

def merger_func(plot_df, cdr_df) -> pd.DataFrame:
    merged_df = pd.merge(plot_df, cdr_df, on='h3_chain', how='left')
    return merged_df

def compute_ctd_vector(sequence: str) -> list:
    ctd_dict = GetProDes(sequence).GetCTD()
    processed_values = []
    for value in ctd_dict.values():
        processed_values.append(float(value))
    return processed_values