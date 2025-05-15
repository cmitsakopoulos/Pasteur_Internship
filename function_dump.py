import os
import sys
import csv
import json
import re
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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
            if headers != expected:
                print("NORMAL csv type, reverting to standard CSV read.", file=sys.stderr)
                return pd.read_csv(filepath) #If its not idiosyncratic
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
    return pd.DataFrame(records) #self explanatory

def parser(df):
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
                antigen_row['antigen_id'] = f"{basic_dict.get('pdb_id')}_antigen"
            #Append the rows to the records lists separately
            cdr_records.append(cdr_row)
            antigen_records.append(antigen_row)
    else: #MODULATE THIS FOR NEW CSVS, CANNOT BE BOTHERED TO ACCOMODATE EVERY USE CASE; MAYBE IN THE FUTURE
            basic_records = []
            for idx, new_row in df.iterrows():
                row = {}
                row['pdb_id'] = new_row['pdb_name']
                row["antigen_seq"] = new_row["antigen_sequence"]
                row["h3_chain"] = new_row["heavy_cdr3"]
                row['l3_chain'] = new_row["light_cdr3"]
                row['database_origin'] = new_row["dataset"]
                basic_records.append(row)
    #Convert the lists to separate dataframes
    antigen_df = pd.DataFrame(antigen_records)
    cdr_df = pd.DataFrame(cdr_records)
    #Return based on use case
    if antigen_df and cdr_df:
        return antigen_df, cdr_df
    else:
        return pd.DataFrame(basic_records)

def calculate_cdr_chars(df):
    #Check for pre-existing columns
    for col in ['h3_gravy', 'h3_pI', 'h3_net_charge_inflamed', 'h3_net_charge_normal', 'l3_gravy', 'l3_pI', 'l3_net_charge_inflamed',  'l3_net_charge_normal']:
        df[col] = float('nan')

    for idx, row in df.iterrows():
        # Heavy chain CDR3
        h3_seq = row.get('h3_chain')
        if not h3_seq:
            print(f"Row {idx} has no H3 sequence.")
        else:
            is_incomplete_3 = 'X' in h3_seq
            clean_seq_h3     = h3_seq.replace('X', '')
            try:
                agn_3 = ProteinAnalysis(clean_seq_h3)
            except Exception as e:
                print(f"Row {idx}: failed to analyze antigen ({e})")
            else:
                val = agn_3.gravy()
                df.at[idx, 'h3_gravy'] = float(f"{val:.5g}")

                val = agn_3.isoelectric_point()
                df.at[idx, 'h3_pI'] = float(f"{val:.5g}")

                val = agn_3.charge_at_pH(7.35)
                df.at[idx, 'h3_net_charge_normal'] = float(f"{val:.5g}")

                val = agn_3.charge_at_pH(5.5)
                df.at[idx, 'h3_net_charge_inflamed'] = float(f"{val:.5g}")
                df.at[idx, 'h3_is_incomplete']       = is_incomplete_3
                df.at[idx, "h3_N_gylcosylation_sites"] = find_N_glycosilation(h3_seq)
                df.at[idx, "h3_O_gylcosylation_sites"] = len(re.findall(r'(?=(?:[ST]{3}))(?!P)', h3_seq))

        # Light chain CDR3
        l3_seq = row.get('l3_chain')
        if not l3_seq:
            print(f"Row {idx} has no L3 sequence.")
        else:
            is_incomplete_2 = 'X' in l3_seq
            clean_seq_l3     = l3_seq.replace('X', '')
            try:
                agn_2 = ProteinAnalysis(clean_seq_l3)
            except Exception as e:
                print(f"Row {idx}: failed to analyze light chain ({e})")
            else:
                val = agn_2.gravy()
                df.at[idx, 'l3_gravy'] = float(f"{val:.5g}")

                val = agn_2.isoelectric_point()
                df.at[idx, 'l3_pI'] = float(f"{val:.5g}")

                val = agn_2.charge_at_pH(7.35)
                df.at[idx, 'l3_net_charge_normal'] = float(f"{val:.5g}")

                val = agn_2.charge_at_pH(5.5)
                df.at[idx, 'l3_net_charge_inflamed'] = float(f"{val:.5g}")
                df.at[idx, 'l3_is_incomplete']       = is_incomplete_2
                df.at[idx, "l3_N_gylcosylation_sites"] = find_N_glycosilation(l3_seq)
                df.at[idx, "l3_O_gylcosylation_sites"] = len(re.findall(r'(?=(?:[ST]{3}))(?!P)', h3_seq))

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

def calculate_antigen_chars(df):
    # If columns exist, continue, else we make them from sratch
    for col in ["antigen_gravy", "antigen_pI", "antigen_net_charge_inflamed", "antigen_net_charge_normal", "antigen_is_incomplete"]:
        df[col] = float('nan')

    for idx, row in df.iterrows():
        print(f"Processing row {idx} sequence data")
        antigen = row.get('antigen_seq')
        if not antigen:
            print(f"Row {idx} has no antigen sequence.")
        else:
            is_incomplete = 'X' in antigen
            clean_seq     = antigen.replace('X', '')
            try:
                agn = ProteinAnalysis(clean_seq)
            except Exception as e:
                print(f"Row {idx}: failed to analyze antigen ({e})")
            else:
                val = agn.gravy()
                df.at[idx, 'antigen_gravy'] = float(f"{val:.5g}")

                val = agn.isoelectric_point()
                df.at[idx, 'antigen_pI'] = float(f"{val:.5g}")

                val = agn.charge_at_pH(7.35)
                df.at[idx, 'antigen_net_charge_normal'] = float(f"{val:.5g}")

                val = agn.charge_at_pH(5.5)
                df.at[idx, 'antigen_net_charge_inflamed'] = float(f"{val:.5g}")
                df.at[idx, 'antigen_is_incomplete']       = is_incomplete

    return df

def inspect_verbose(df):
    for idx, row in df.iterrows():
        print(f"\nRow {idx}:")
        for key, val in row.items():
            print(f"  {key}: {val}")

def inspect_summary(df):
    print(f"Provided DataFrame contains {len(df)} records")