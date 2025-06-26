import os
import sys
import csv
import json
import re
import shutil
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
from propy import Autocorrelation

def prepare_charges_hydrophobicity(): #Better to keep things fresh at each run, rather than rely on a predefined dictionary
    aas = "ACDEFGHIKLMNPQRSTVWY"
    aa_list = [i for i in aas]
    net_charge_at_pH735 = {}
    net_charge_at_pH550 = {}
    hydrophobicity_scale = {}
    for index in aa_list:
        aa_prep = ProteinAnalysis(index)
        net_charge_at_pH735[index] = aa_prep.charge_at_pH(7.35)
        net_charge_at_pH550[index] = aa_prep.charge_at_pH(5.5)
        hydrophobicity_scale[index] = ProtParamData.kd[index]
    return net_charge_at_pH550, net_charge_at_pH735, hydrophobicity_scale

def calc_autocorrelations_mini(seq, blood_geary_charge: dict, inflamed_geary_charge: dict, hydrophobicity_scale: dict):
    #Remember the precalculated dicts above (def prepare_charges_hydrophobicity)?
    ph735 = Autocorrelation.CalculateEachGearyAuto(
        seq,
        blood_geary_charge,
        "net_charges",
    )
    ph550 = Autocorrelation.CalculateEachGearyAuto(
        seq,
        inflamed_geary_charge,
        "net_charges"
    )
    hydrophobicity = Autocorrelation.CalculateEachGearyAuto(
        seq,
        hydrophobicity_scale,
        "hydrophobicity"
    )
    return list(ph550.values()), list(ph735.values()), list(hydrophobicity.values())
#For visual purposes, Ive split the functions into a separate file so that it is easier to debug...only personal, literally-obviously zero difference...
csv.field_size_limit(sys.maxsize)
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

"""
Introduced new function to check for the dataframe "health", if it is from the pipeline, 
then continue correctly.
"""
def parser(df):
    parsed_cdr_cols = {'h3_chain', 'l3_chain', 'pdb_id'}
    parsed_antigen_cols = {'antigen_seq', 'database_origin', "last_update", "heavy_host_organism_name"}
    if parsed_cdr_cols.issubset(df.columns):
        print("Skipping nascent parsing, DataFrame already contains CDR columns.")
        empty_antigen_df = pd.DataFrame(...)
        return empty_antigen_df, df

    if parsed_antigen_cols.issubset(df.columns):
        print("Skipping nascent parsing, DataFrame already contains Antigen columns.")
        empty_cdr_df = pd.DataFrame(...)
        return df, empty_cdr_df
    
    cdr_records = []
    antigen_records = []
    if 'basic' in df.columns: # IF IDIOSYNCRATIC
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
                # Split the combined l3_h3 string into separate chains
                l3_h3 = basic_dict.get('l3_h3', '')
                parts = l3_h3.split('_', 1) #max one split
                cdr_row['l3_chain'] = parts[0] if len(parts) > 0 else None
                cdr_row['h3_chain'] = parts[1] if len(parts) > 1 else None
                antigen_row['corresponding_pdb_antibody'] = f"{basic_dict.get('pdb_id')}"
            cdr_records.append(cdr_row)
            antigen_records.append(antigen_row)
    else: #MODULATE THIS FOR NEW CSVS, CANNOT BE BOTHERED TO ACCOMODATE EVERY USE CASE; MAYBE IN THE FUTURE
            cdr_records = []
            antigen_records = []
            for idx, new_row in df.iterrows():
                cdr_row = {}
                antigen_row = {}
                antigen_row['corresponding_pdb_antibody'] = new_row['pdb_name']
                cdr_row['pdb_id'] = new_row['pdb_name']
                antigen_row["antigen_seq"] = new_row["antigen_sequence"]
                cdr_row["h3_chain"] = new_row["heavy_cdr3"]
                cdr_row['l3_chain'] = new_row["light_cdr3"]
                cdr_row['database_origin'] = new_row["dataset"]
                antigen_row['database_origin'] = new_row["dataset"]
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

def calculate_cdr_chars(df):
    for column in ["h3_geary_hydrophobicity", "h3_blood_geary_charge", "h3_inflamed_geary_charge", "l3_blood_geary_charge", "l3_inflamed_geary_charge", "l3_geary_hydrophobicity"]:
        df[column] = None
    dict_inflamed, dict_normal, dict_hydro = prepare_charges_hydrophobicity()
    for idx, row in df.iterrows():
        # Heavy chain CDR3
        h3_seq = row.get('h3_chain')
        if not h3_seq:
            print(f"Row {idx} has no H3 sequence.")
            return None
        else:
            is_incomplete_3 = int('X' in h3_seq)
            clean_seq_h3     = h3_seq.replace('X', '')
            try:
                agn_3 = ProteinAnalysis(clean_seq_h3)
            except Exception as e:
                print(f"Row {idx}: failed to analyze antigen ({e})")
            else:
                inflamed_h3, normal_h3, hydrophobicity_h3 = calc_autocorrelations_mini(clean_seq_h3, dict_normal, dict_inflamed, dict_hydro)
                df.at[idx, 'h3_geary_hydrophobicity'] = hydrophobicity_h3
                val = agn_3.isoelectric_point()
                df.at[idx, 'h3_pi'] = float(f"{val:.5g}")
                df.at[idx, 'h3_blood_geary_charge'] = normal_h3
                df.at[idx, 'h3_inflamed_geary_charge'] = inflamed_h3
                df.at[idx, 'h3_is_incomplete']       = is_incomplete_3
                df.at[idx, "h3_n_glycosilation_sites"] = find_N_glycosilation(h3_seq)
                df.at[idx, "h3_o_glycosilation_sites"] = len(re.findall(r'(?=(?:[ST]{3}))(?!P)', h3_seq))
                df.at[idx, "h3_chain"] = clean_seq_h3

        # Light chain CDR3
        l3_seq = row.get('l3_chain')
        if not l3_seq:
            print(f"Row {idx} has no L3 sequence.")
        else:
            is_incomplete_2 = int('X' in l3_seq)
            clean_seq_l3     = l3_seq.replace('X', '')
            try:
                agn_2 = ProteinAnalysis(clean_seq_l3)
            except Exception as e:
                print(f"Row {idx}: failed to analyze light chain ({e})")
            else:
                inflamed_l3, normal_l3, hydrophobicity_l3 = calc_autocorrelations_mini(clean_seq_l3, dict_normal, dict_inflamed, dict_hydro)
                df.at[idx, 'l3_geary_hydrophobicity'] = hydrophobicity_l3

                val = agn_2.isoelectric_point()
                df.at[idx, 'l3_pi'] = float(f"{val:.5g}")
                df.at[idx, 'l3_blood_geary_charge'] = normal_l3
                df.at[idx, 'l3_inflamed_geary_charge'] = inflamed_l3
                df.at[idx, 'l3_is_incomplete']       = is_incomplete_2
                df.at[idx, "l3_n_glycosilation_sites"] = find_N_glycosilation(l3_seq)
                df.at[idx, "l3_o_glycosilation_sites"] = len(re.findall(r'(?=(?:[ST]{3}))(?!P)', h3_seq))
                df.at[idx, "l3_chain"] = clean_seq_l3
    return df

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

def calculate_antigen_chars(df):#Copy of the antibody related chemical chars function
    for column in  ["antigen_geary_hydrophobicity", "antigen_blood_geary_charge", "antigen_inflamed_geary_charge"]:
        df[column] = None
    dict_inflamed, dict_normal, dict_hydro = prepare_charges_hydrophobicity() #Why not call it in the Step class which handles this function? Antigen/CDR are handled separate, therefore no difference in running it here or the separate Step(s)
    for idx, row in df.iterrows():
        antigen = row.get('antigen_seq')
        if not antigen:
            print(f"Row {idx} has no antigen sequence.")
            return None
        else:
            is_incomplete = int('X' in antigen)
            clean_seq     = antigen.replace('X', '')
            try:
                agn = ProteinAnalysis(clean_seq)
            except Exception as e:
                print(f"Row {idx}: failed to analyze antigen ({e})")
            else:
                inflamed_ant, normal_ant, hydrophobicity_ant = calc_autocorrelations_mini(clean_seq, dict_normal, dict_inflamed, dict_hydro)
                df.at[idx, 'antigen_geary_hydrophobicity'] = hydrophobicity_ant
                val = agn.isoelectric_point()
                df.at[idx, 'antigen_pi'] = float(f"{val:.5g}")
                df.at[idx, 'antigen_blood_geary_charge'] = normal_ant
                df.at[idx, 'antigen_inflamed_geary_charge'] = inflamed_ant
                df.at[idx, 'antigen_is_incomplete']       = is_incomplete
                df.at[idx, "antigen_seq"] = clean_seq
    return df

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

def clear_dir(path: str):
    if not os.path.isdir(path):
        return
    for name in os.listdir(path):
        full = os.path.join(path, name)
        if os.path.isfile(full) or os.path.islink(full):
            os.remove(full)
        elif os.path.isdir(full):
            shutil.rmtree(full)