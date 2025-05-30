from sqlalchemy import create_engine
import pandas as pd

sql_engine = create_engine("postgresql://chrismitsacopoulos:password@localhost/pasteurdb")
df = pd.read_sql_query("SELECT * FROM training_dataset", sql_engine)

# Identify duplicate sequence entries based on h3_chain, l3_chain, and antigen_sequence
dup_mask = df.duplicated(subset=['h3_chain', 'l3_chain', 'antigen_seq'], keep=False)
duplicates = df[dup_mask]
print(f"Number of duplicate rows: {duplicates.shape[0]}")
if not duplicates.empty:
    print("Duplicate entries found:")
    print(duplicates[['h3_chain', 'l3_chain', 'antigen_seq']])

