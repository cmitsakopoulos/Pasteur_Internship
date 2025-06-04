from sqlalchemy import create_engine
import pandas as pd
import numpy as np
import time
from sklearn.preprocessing import OneHotEncoder, FunctionTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline as SeqPipeline
from sklearn.model_selection import KFold
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping


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

# Standardize numeric features for neural network
x_scaler = StandardScaler()
X_numeric =x_scaler.fit_transform(X_numeric)


Y_numeric = df[["antigen_isoelectric", "antigen_gravy", "antigen_net_charge_normal", "antigen_net_charge_inflamed"]]
y_scaler = StandardScaler()
Y_numeric = y_scaler.fit_transform(Y_numeric)

# Encode L3 chain sequences and combine with H3
X_l3 = antibody_encoder.transform(df[["l3_chain"]])
X_seq = np.hstack([X_h3, X_l3])

# Combine sequence encodings and numeric features
X = np.hstack([X_seq, X_numeric])

# --- Build and evaluate with Keras ---
def build_model(input_dim, output_dim):
    model = Sequential([
        Dense(100, activation='relu', input_shape=(input_dim,)),
        Dense(100, activation='relu'),
        Dense(output_dim, activation='linear'),
    ])
    model.compile(optimizer=Adam(), loss='mse', metrics=['mae'])
    return model

# K-Fold cross-validation
kf = KFold(n_splits=5, shuffle=True, random_state=69)
all_histories = []
start_time = time.time()
for fold, (train_idx, val_idx) in enumerate(kf.split(X), 1):
    print(f"\n--- Fold {fold} ---")
    X_train, X_val = X[train_idx], X[val_idx]
    Y_train, Y_val = Y_numeric[train_idx], Y_numeric[val_idx]
    
    model = build_model(X.shape[1], Y_numeric.shape[1])
    history = model.fit(
        X_train, Y_train,
        validation_data=(X_val, Y_val),
        epochs=50,
        batch_size=32,
        callbacks=[EarlyStopping(patience=5, restore_best_weights=True)],
        verbose=1
    )
    val_loss, val_mae = model.evaluate(X_val, Y_val, verbose=0)
    print(f"Validation loss: {val_loss:.4f}, MAE: {val_mae:.4f}")
    all_histories.append(history)
end_time = time.time()
print(f"\nTotal training time: {end_time - start_time:.2f} seconds")
