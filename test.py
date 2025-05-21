import pandas as pd
import re

antigens = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/antigen_20250521_182332.csv")
cdrs = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/cdr_20250521_182332.csv")

def find_cdr_ids(txt):
    txt = str(txt)
    mask = cdrs["pdb_id"].str.contains(txt, na=False)
    return cdrs.loc[mask, "cdr_computed_id"].tolist()

antigens["matched_cdr_ids"] = (
    antigens["corresponding_pdb_antibody"]
      .map(find_cdr_ids)
)

print(antigens["matched_cdr_ids"])