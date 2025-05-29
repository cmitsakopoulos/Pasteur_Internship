import pandas as pd
from typing import Dict

# Utility to convert PDB fields into lists
def to_list(x):
    if pd.isna(x):
        return []
    if isinstance(x, str):
        return [item.strip() for item in x.split(',') if item.strip()]
    if isinstance(x, list):
        return x
    return []

class ComputeRelationships:
    def process(self, data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        new_data: Dict[str, pd.DataFrame] = {}
        antigens = data["antigen"].copy()
        cdrs = data["cdr"].copy()

        # Convert PDB columns to lists
        antigens["pdb_list"] = antigens["corresponding_pdb_antibody"].apply(to_list)
        cdrs["pdb_list"]     = cdrs["pdb_id"].apply(to_list)

        # Explode to one row per pdb entry
        ag_exp  = antigens.explode("pdb_list")[["antigen_computed_id", "pdb_list"]]
        cdr_exp = cdrs.explode("pdb_list")[["cdr_computed_id",     "pdb_list"]]

        # Merge on pdb_list and dedupe
        pairs = (
            pd.merge(ag_exp, cdr_exp, on="pdb_list", how="inner")
              [["antigen_computed_id", "cdr_computed_id"]]
              .drop_duplicates()
              .reset_index(drop=True)
        )

        # Assign outputs, dropping helper column
        new_data["antigen"]       = antigens.drop(columns=["pdb_list"])
        new_data["cdr"]           = cdrs.drop(columns=["pdb_list"])
        new_data["relationships"] = pairs
        return new_data