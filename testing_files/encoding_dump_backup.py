from sqlalchemy import create_engine
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor, StackingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import time
from sklearn.preprocessing import OneHotEncoder, FunctionTransformer
from sklearn.pipeline import Pipeline as SeqPipeline


sql_engine = create_engine("postgresql://chrismitsacopoulos:password@localhost/pasteurdb")
df = pd.read_sql_query("SELECT * FROM training_dataset", sql_engine)

"""
Big problem-> For one hot encoding to actually work, we need to make all sequences in one column be the same length. Therefore, the split_pad function introduces X characters which are NOT in the alphabet (see amino_acids) and OneHotEncoder has been instructed to IGNORE them when creating categories and mapping out sequence features...
"""
amino_acids = list("ACDEFGHIKLMNPQRSTVWY") #Ported from the pipeline...hopefully no mistakes happened there or this is completely wrong!!!
max_len = df[["h3_chain", "l3_chain"]].applymap(len).max().max() #Use pandas to identify the longest possible entry for both columns, such that when performing one hot encoding, categories will be built for the max length sequence, smaller sequences are handled appropriately by onehotencoding from sklearn...I hope :)

antigen_max_len = df["antigen_sequence"].str.len().max()

# generic split-and-pad factory
def split_pad_factory(length):
    def split_pad(seqs):
        return [list(seq.ljust(length, "X")[:length]) for seq in seqs.ravel()]
    return split_pad

#Sklearn gives us the ability to make a pipeline to do the encoding steps for us, this is it for X sequence variables...
#Do the one hot encoding separately, such that based on sequence length we create categories that are length specific, to avoid padding. Considering that both antigens and antiboides share the same amion acid basis, categories should be created with a common "logic base" per say.
encoded_sequences = SeqPipeline([
    ("split", FunctionTransformer(split_pad_factory(max_len), validate=False)),
    ("onehot", OneHotEncoder(
        categories=[amino_acids] * max_len,
        handle_unknown="ignore", #This will ignore the added/padding X characters
        sparse_output=False
    ))
])

#Re introduce the encoded sequences back into the dataframe
encoded_sequences.fit(df[["h3_chain", "l3_chain"]])
X_h3 = encoded_sequences.transform(df[["h3_chain"]])
X_l3 = encoded_sequences.transform(df[["l3_chain"]])
X_seq = np.hstack([X_h3, X_l3])

antigen_encoder = SeqPipeline([
    ("split", FunctionTransformer(split_pad_factory(antigen_max_len), validate=False)),
    ("onehot", OneHotEncoder(
        categories=[amino_acids] * antigen_max_len,
        handle_unknown="ignore",
        sparse_output=False
    ))
])

antigen_encoder.fit(df[["antigen_sequence"]])
Y_antigen = antigen_encoder.transform(df[["antigen_sequence"]])

# Unpack list-valued "geary" columns into 2D arrays
h3_geary = np.vstack(df['h3_geary_hydrophobicity'].values)
l3_geary = np.vstack(df['l3_geary_hydrophobicity'].values)
h3_blood_charge = np.vstack(df['h3_blood_geary_charge'].values)
l3_blood_charge = np.vstack(df['l3_blood_geary_charge'].values)

# Scalar numeric columns
scalar_cols = [
    "h3_pi", "l3_pi",
    "h3_n_glycosilation_sites", "h3_o_glycosilation_sites",
    "l3_n_glycosilation_sites", "l3_o_glycosilation_sites"
]
scalar_numeric = df[scalar_cols].values

# Combine all numeric parts into X_numeric
X_numeric = np.hstack([h3_geary, l3_geary, h3_blood_charge, l3_blood_charge, scalar_numeric])

Y_numeric = df[["antigen_isoelectric", "antigen_geary_hydrophobicity", "antigen_blood_geary_charge", "antigen_blood_geary_charge"]]

X = np.hstack([X_numeric, X_seq])  
Y = np.hstack([Y_antigen, Y_numeric])