import os
import sys
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from propy import CTD
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

from app_components.config import COMPUTATION_DIR

def compute_ctd_vector(sequence: str) -> np.ndarray:
    if not isinstance(sequence, str) or len(sequence) < 2:
        return np.array([])
    try:
        c_features = CTD.CalculateC(sequence)
        t_features = CTD.CalculateT(sequence)
        d_features = CTD.CalculateD(sequence)
        all_features = {**c_features, **t_features, **d_features}
        return np.array(list(all_features.values()))
    except Exception as e:
        print(f"Error computing CTD for sequence '{sequence}': {e}")
        return np.array([])

def calculate_blosum_score(seq1: str, seq2:str) -> float:
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM45")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -1.0
    try:
        score = aligner.score(seq1, seq2)
        return score
    except Exception:
        return 0.0

df_main = pd.read_csv(os.path.join(COMPUTATION_DIR, "cdr.csv"))
df_main.dropna(subset=['h3_chain'], inplace=True)
df_main = df_main.head(100).reset_index(drop=True)

print("-> Calculating CTD vectors for each sequence...")
df_main['ctd_vector'] = df_main['h3_chain'].apply(compute_ctd_vector)

df_filtered = df_main[df_main['ctd_vector'].apply(len) > 0].reset_index(drop=True)

if df_filtered.empty:
    print("\nError: No valid CTD vectors were generated. Halting execution.")
    print("Please check the output above for errors from the 'compute_ctd_vector' function.")
    sys.exit(1)

ctd_feature_matrix = np.vstack(df_filtered['ctd_vector'].values)

print("-> Computing Manhattan distance matrix from CTD vectors...")
ctd_distance_matrix = pairwise_distances(ctd_feature_matrix, metric='manhattan')

print("-> Calculating pairwise BLOSUM similarity matrix...")
num_sequences = len(df_filtered)
sequences = df_filtered['h3_chain'].tolist()
blosum_similarity_matrix = np.zeros((num_sequences, num_sequences))

for i in range(num_sequences):
    for j in range(i, num_sequences):
        score = calculate_blosum_score(sequences[i], sequences[j])
        blosum_similarity_matrix[i, j] = blosum_similarity_matrix[j, i] = score

print("-> Converting BLOSUM similarity to distance and normalizing with MinMaxScaler...")
blosum_distance_matrix = np.max(blosum_similarity_matrix) - blosum_similarity_matrix

scaler = MinMaxScaler()
normalized_blosum_distance_matrix = scaler.fit_transform(blosum_distance_matrix)


print("\n--- Computation Complete ---")
print(f"Shape of CTD Manhattan Distance Matrix: {ctd_distance_matrix.shape}")
print(f"Shape of Normalized BLOSUM Distance Matrix: {normalized_blosum_distance_matrix.shape}")
print(f"Data filtered to {num_sequences} sequences with valid features.")