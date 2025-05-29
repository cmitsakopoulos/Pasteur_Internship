from sqlalchemy import create_engine
import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder, FunctionTransformer
from sklearn.pipeline import Pipeline as SeqPipeline


sql_engine = create_engine("postgresql://chrismitsacopoulos:password@localhost/pasteurdb")
df = pd.read_sql_query("SELECT * FROM training_dataset", sql_engine)

"""
Big problem-> For one hot encoding to actually work, we need to make all sequences in one column be the same length. Therefore, the split_pad function introduces X characters which are NOT in the alphabet (see amino_acids) and OneHotEncoder has been instructed to IGNORE them when creating categories and mapping out sequence features...
"""
amino_acids = list("ACDEFGHIKLMNPQRSTVWY") #Ported from the pipeline...hopefully no mistakes happened there or this is completely wrong!!!
# compute max sequence length for both antibody chains using Series.map()
max_len = max(df["h3_chain"].map(len).max(), df["l3_chain"].map(len).max())

antigen_max_len = df["antigen_sequence"].str.len().max()

# generic split-and-pad factory
def split_pad_factory(length):
    def split_pad(seqs):
        seqs_arr = np.array(seqs).ravel()
        return [list(seq.ljust(length, "X")[:length]) for seq in seqs_arr]
    return split_pad

#Sklearn gives us the ability to make a pipeline to do the encoding steps for us, this is it for X sequence variables...

antibody_encoder = SeqPipeline([
    ("split", FunctionTransformer(split_pad_factory(max_len), validate=False)),
    ("onehot", OneHotEncoder(
        categories=[amino_acids] * max_len,
        handle_unknown="ignore", #This will ignore the added/padding X characters
        sparse_output=False
    ))
])

#Re introduce the encoded sequences back into the dataframe
antibody_encoder.fit(df[["h3_chain", "l3_chain"]])
X_h3 = antibody_encoder.transform(df[["h3_chain"]])
#X_l3 = antibody_encoder.transform(df[["l3_chain"]])
#X_seq = np.hstack([X_h3, X_l3])

antigen_encoder = SeqPipeline([
    ("split", FunctionTransformer(split_pad_factory(antigen_max_len), validate=False)),
    ("onehot", OneHotEncoder(
        categories=[amino_acids] * antigen_max_len,
        handle_unknown="ignore",
        sparse_output=False
    ))
])

antigen_encoder.fit(df[["antigen_sequence"]])
Y_seq = antigen_encoder.transform(df[["antigen_sequence"]])

#    "l3_n_gylcosylation_sites", "l3_o_gylcosylation_sites",
#   "l3_net_charge_inflamed", "l3_net_charge_normal",
#  "l3_gravy", "l3_pi",
X_numeric = df[[
    "h3_gravy", "h3_pi",
    "h3_n_gylcosylation_sites", "h3_o_gylcosylation_sites",
    "h3_net_charge_inflamed", "h3_net_charge_normal"
]].values

Y_numeric = df[["antigen_isoelectric", "antigen_gravy", "antigen_net_charge_normal", "antigen_net_charge_inflamed"]]
