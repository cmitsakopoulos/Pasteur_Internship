from sqlalchemy import create_engine
import pandas as pd
import numpy as np
import time
from sklearn.preprocessing import OneHotEncoder, FunctionTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline as SeqPipeline
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping

sql_engine = create_engine("postgresql://chrismitsacopoulos:password@localhost/pasteurdb")
df = pd.read_sql_query("SELECT * FROM training_dataset", sql_engine)

"""
Big problem-> For one hot encoding to actually work, we need to make all sequences in one column be the same length. Therefore, the split_pad function introduces X characters which are NOT in the alphabet (see amino_acids) and OneHotEncoder has been instructed to IGNORE them when creating categories and mapping out sequence features...
"""
amino_acids = list("ACDEFGHIKLMNPQRSTVWY") ax() 

antigen_max_len = df["antigen_sequence"].str.len().max()

h3_geary = np.vstack(df['h3_geary_hydrophobicity'].values)
l3_geary = np.vstack(df['l3_geary_hydrophobicity'].values)
h3_blood_charge = np.vstack(df['h3_blood_geary_charge'].values)
l3_blood_charge = np.vstack(df['l3_blood_geary_charge'].values)
h3_inflamed_charge = np.vstack(df['h3_inflamed_geary_charge'].values)
l3_inflamed_charge = np.vstack(df['l3_inflamed_geary_charge'].values)

h3_mask = np.any(h3_geary != 0.0, axis=0)
h3_geary = h3_geary[:, h3_mask]

l3_mask = np.any(l3_geary != 0.0, axis=0)
l3_geary = l3_geary[:, l3_mask]

h3_blood_mask = np.any(h3_blood_charge != 0.0, axis=0)
h3_blood_charge = h3_blood_charge[:, h3_blood_mask]

l3_blood_mask = np.any(l3_blood_charge != 0.0, axis=0)
l3_blood_charge = l3_blood_charge[:, l3_blood_mask]

h3_inflamed_mask = np.any(h3_inflamed_charge != 0.0, axis=0)
h3_inflamed_charge = h3_inflamed_charge[:, h3_inflamed_mask]

l3_inflamed_mask = np.any(l3_inflamed_charge != 0.0, axis=0)
l3_inflamed_charge = l3_inflamed_charge[:, l3_inflamed_mask]

scalar_cols = [
    "h3_pi", "l3_pi",
    "h3_n_glycosilation_sites", "h3_o_glycosilation_sites",
    "l3_n_glycosilation_sites", "l3_o_glycosilation_sites"
]
scalar_numeric = df[scalar_cols].values

X_numeric = np.hstack([
    h3_geary, l3_geary,
    h3_blood_charge, l3_blood_charge,
    h3_inflamed_charge, l3_inflamed_charge,
    scalar_numeric
])

antigen_geary = np.vstack(df['antigen_geary_hydrophobicity'].values)
antigen_blood_charge = np.vstack(df['antigen_blood_geary_charge'].values)
antigen_inflamed_charge = np.vstack(df['antigen_inflamed_geary_charge'].values)

antigen_geary_mask = np.any(antigen_geary != 0.0, axis=0)
antigen_geary = antigen_geary[:, antigen_geary_mask]

antigen_blood_mask = np.any(antigen_blood_charge != 0.0, axis=0)
antigen_blood_charge = antigen_blood_charge[:, antigen_blood_mask]

antigen_inflamed_mask = np.any(antigen_inflamed_charge != 0.0, axis=0)
antigen_inflamed_charge = antigen_inflamed_charge[:, antigen_inflamed_mask]

antigen_iso = df['antigen_isoelectric'].values.reshape(-1, 1)

Y_num_array = np.hstack([antigen_iso, antigen_geary, antigen_blood_charge, antigen_inflamed_charge]).astype(np.float32)

X = X_numeric  
X = X.astype(np.float32)
Y = Y_num_array  

def build_model(input_dim, output_dim):
    model = Sequential([
        Dense(100, activation='relu', input_shape=(input_dim,)),
        Dense(100, activation='relu'),
        Dense(output_dim, activation='linear'),
    ])
    model.compile(optimizer=Adam(), loss='mse', metrics=['mae'])
    return model

kf = KFold(n_splits=5, shuffle=True, random_state=69)
all_histories = []
start_time = time.time()
for fold, (train_idx, val_idx) in enumerate(kf.split(X), 1):
    print(f"\n--- Fold {fold} ---")
    X_train, X_val = X[train_idx], X[val_idx]
    Y_train, Y_val = Y_num_array[train_idx], Y_num_array[val_idx]

    model = build_model(X.shape[1], Y_num_array.shape[1])

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
    y_pred = model.predict(X_val)
    r2 = r2_score(Y_val, y_pred)
    print(f"R² on validation numeric targets: {r2:.4f}")
    all_histories.append(history)
end_time = time.time()
print(f"\nTotal training time: {end_time - start_time:.2f} seconds")