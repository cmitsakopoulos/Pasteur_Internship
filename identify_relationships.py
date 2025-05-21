import pandas as pd

antigens = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/antigen_20250521_182332.csv")
cdrs = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/cdr_20250521_182332.csv")

def normalise(col):
    return (
        col.apply(
            lambda x: [] if pd.isna(x)
            else x if isinstance(x, list)
            else [x]                     
        )
    )

a_long = (antigens
          .assign(corresponding_pdb_antibody=normalise(antigens["corresponding_pdb_antibody"]))
          .explode("corresponding_pdb_antibody"))        

c_long = (cdrs
          .assign(pdb_id=normalise(cdrs["pdb_id"]))
          .explode("pdb_id"))

pairs = (a_long
         .merge(c_long,
                left_on="corresponding_pdb_antibody",
                right_on="pdb_id",
                how="inner")[["antigen_computed_id", "cdr_computed_id"]]
         .drop_duplicates())

print(pairs)