import os
import sys
import pandas as pd
from propy.PyPro import GetProDes


project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

from app_components.config import COMPUTATION_DIR

df_main = pd.read_csv(os.path.join(COMPUTATION_DIR, "cdr.csv"))
df_main.dropna(subset=['h3_chain'], inplace=True)
df_main = df_main.head(100).reset_index(drop=True)

def compute_ctd_vector(sequence: str) -> list:
    ctd_dict = GetProDes(sequence).GetCTD()
    processed_values = []
    for value in ctd_dict.values():
        if isinstance(value, (list, tuple)):
            processed_values.extend([float(item) for item in value])
        else:
            processed_values.append(float(value))
    return processed_values

df_main['h3_ctd'] = df_main['h3_chain'].apply(compute_ctd_vector)

print(df_main["h3_ctd"])