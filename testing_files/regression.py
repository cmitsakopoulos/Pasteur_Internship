from sqlalchemy import create_engine
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor, StackingRegressor
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, cross_val_score, KFold
from sklearn.model_selection import cross_validate
import time
from sklearn.preprocessing import OneHotEncoder, FunctionTransformer
from sklearn.pipeline import Pipeline as SeqPipeline
from xgboost import XGBRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.multioutput import MultiOutputRegressor


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

#    "l3_n_glycosilation_sites", "l3_o_glycosilation_sites",
#   "l3_net_charge_inflamed", "l3_net_charge_normal",
#  "l3_gravy", "l3_pi",
X_numeric = df[[
    "h3_gravy", "h3_pi",
    "h3_n_glycosilation_sites", "h3_o_glycosilation_sites",
    "h3_net_charge_inflamed", "h3_net_charge_normal"
]].values

x_scaler = StandardScaler()
X_numeric =x_scaler.fit_transform(X_numeric)

Y_numeric = df[["antigen_isoelectric", "antigen_gravy", "antigen_net_charge_normal", "antigen_net_charge_inflamed"]]
y_scaler = StandardScaler()
Y_numeric = y_scaler.fit_transform(Y_numeric)

stack_estimators = [
    ('cont_xgb', Pipeline([
        ('scale', StandardScaler()),
        ('hgb', MultiOutputRegressor(
            HistGradientBoostingRegressor(
                max_iter=1000,
                learning_rate=0.05,
                max_depth=10,
                l2_regularization=1.0,
                random_state=42
            )
        ))
    ])),
    ('seq_ridge', RidgeCV(
        alphas=[0.01, 0.1, 1.0, 10.0, 100.0],
        scoring='neg_mean_squared_error',
        cv=5
    ))
]

stack = StackingRegressor(
    estimators=stack_estimators,
    final_estimator=Ridge(),   
    passthrough=True         
)



# Extract separate modality pipelines
cont_xgb_pipeline = dict(stack_estimators)['cont_xgb']
seq_ridge_pipeline = dict(stack_estimators)['seq_ridge']


# -- Cross-validation evaluation (multiple metrics at once) --
kf = KFold(n_splits=5, shuffle=True, random_state=42)
scoring = {
    'rmse': 'neg_root_mean_squared_error',
    'r2': 'r2'
}

# Numeric features CV
start_num = time.time()
cv_num = cross_validate(
    cont_xgb_pipeline,
    X_numeric, Y_numeric,
    cv=kf,
    scoring=scoring,
    return_train_score=False
)
elapsed_num = time.time() - start_num
rmse_num = np.sqrt(-cv_num['test_rmse'])
r2_num = cv_num['test_r2']
print(f"Numeric CV RMSE: {rmse_num.mean():.3f} ± {rmse_num.std():.3f}")
print(f"Numeric CV R²:  {r2_num.mean():.3f} ± {r2_num.std():.3f}")
print(f"Numeric CV time: {elapsed_num:.2f}s")

# Sequence features CV
start_seq = time.time()
cv_seq = cross_validate(
    seq_ridge_pipeline,
    X_h3, Y_seq,
    cv=kf,
    scoring=scoring,
    return_train_score=False
)
elapsed_seq = time.time() - start_seq
rmse_seq = np.sqrt(-cv_seq['test_rmse'])
r2_seq = cv_seq['test_r2']
print(f"Sequence CV RMSE: {rmse_seq.mean():.3f} ± {rmse_seq.std():.3f}")
print(f"Sequence CV R²:  {r2_seq.mean():.3f} ± {r2_seq.std():.3f}")
print(f"Sequence CV time: {elapsed_seq:.2f}s")