from sqlalchemy import create_engine
import numpy as np
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Flatten, Dense
from tensorflow.keras.utils import to_categorical
import pandas as pd
from sklearn.model_selection import train_test_split
import time


sql_engine = create_engine("postgresql://chrismitsacopoulos:password@localhost/pasteurdb")
df = pd.read_sql_query("SELECT * FROM training_dataset", sql_engine)

AA = list("ACDEFGHIKLMNPQRSTVWY-")
mapping = {aa: idx for idx, aa in enumerate(AA)}
# fixed lengths for H3 and L3 chains
max_h3 = df['h3_chain'].str.len().max()
max_l3 = df['l3_chain'].str.len().max()
seq_len = max_h3 + max_l3

# pad and concatenate antibody sequences
df['concat_pad'] = df['h3_chain'].str.ljust(max_h3, '-') + df['l3_chain'].str.ljust(max_l3, '-')
# convert to index arrays and one-hot encode
seq_indices = np.array([[mapping[aa] for aa in seq] for seq in df['concat_pad']])
X_seq = to_categorical(seq_indices, num_classes=len(AA))

# numeric antibody features
X_num = df[
    [
      "h3_gravy", "l3_gravy", "h3_pi", "l3_pi",
      "h3_n_gylcosylation_sites", "h3_o_gylcosylation_sites",
      "l3_n_gylcosylation_sites", "l3_o_gylcosylation_sites",
      "l3_net_charge_inflamed", "l3_net_charge_normal",
      "h3_net_charge_inflamed", "h3_net_charge_normal"
    ]
].values

# --- One-hot encoding for antigen sequence targets and numeric targets ---
# pad antigen sequences
max_ag = df['antigen_sequence'].str.len().max()
df['ag_pad'] = df['antigen_sequence'].str.ljust(max_ag, '-')
# convert to index arrays and one-hot encode
ag_indices = np.array([[mapping[aa] for aa in seq] for seq in df['ag_pad']])
Y_seq = to_categorical(ag_indices, num_classes=len(AA))
# numeric antigen targets
Y_num = df[
    [
      "antigen_isoelectric", "antigen_gravy",
      "antigen_net_charge_normal", "antigen_net_charge_inflamed"
    ]
].values

# Split into training and hold-out test sets for both inputs and outputs
X_seq_train, X_seq_test, X_num_train, X_num_test, Y_seq_train, Y_seq_test, Y_num_train, Y_num_test = train_test_split(
    X_seq, X_num, Y_seq, Y_num, test_size=0.1, random_state=42
)

# --- Prepare inputs for a basic MLP ---
# Example: flatten sequence one-hot and concatenate with numeric features
X_train = np.concatenate([
    X_seq_train.reshape((X_seq_train.shape[0], -1)),
    X_num_train
], axis=1)
X_test = np.concatenate([
    X_seq_test.reshape((X_seq_test.shape[0], -1)),
    X_num_test
], axis=1)

# For demonstration, we'll train on numeric antigen targets
y_train = Y_num_train
y_test = Y_num_test

# --- Build a simple feedforward network ---
model = Sequential([
    Flatten(input_shape=(X_train.shape[1],)),
    Dense(64, activation='relu'),
    Dense(64, activation='relu'),
    Dense(y_train.shape[1], activation='linear')
])

model.compile(optimizer='adam', loss='mse', metrics=['mae'])

# --- Train and evaluate ---
history = model.fit(
    X_train, y_train,
    validation_data=(X_test, y_test),
    epochs=10,
    batch_size=32
)

loss, mae = model.evaluate(X_test, y_test)
print(f"Test loss: {loss:.4f}, Test MAE: {mae:.4f}")