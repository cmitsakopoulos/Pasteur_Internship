import numpy as np
import pandas as pd
from propy.PyPro import GetProDes
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

df = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/cdr.csv")

def compute_ctd_vector(sequence: str) -> np.ndarray:
    try:
        if not isinstance(sequence, str) or len(sequence) < 2:
            return np.array([])
        descriptor = GetProDes(sequence)
        ctd_dict = descriptor.GetCTD()
        return np.array(list(ctd_dict.values()))
    except Exception:
        return np.array([])

def calculate_blosum_score(seq1: str, seq2: str) -> float:
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM45")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    try:
        best_score = aligner.score(seq1, seq2)
        return best_score
    except Exception:
        return 0.0

df.dropna(subset=['h3_chain'], inplace=True)
df['ctd_vector'] = df['h3_chain'].apply(compute_ctd_vector)
df = df[df['ctd_vector'].apply(len) > 0]
df.reset_index(drop=True, inplace=True)

feature_matrix = np.vstack(df['ctd_vector'].values)

gram_matrix = feature_matrix @ feature_matrix.T

num_sequences = len(df)
blosum_matrix = np.zeros((num_sequences, num_sequences))
sequences = df['h3_chain'].tolist()

for i in range(num_sequences):
    for j in range(i, num_sequences):
        score = calculate_blosum_score(sequences[i], sequences[j])
        blosum_matrix[i, j] = score
        blosum_matrix[j, i] = score

def min_max_scale(matrix):
    min_val = np.min(matrix)
    max_val = np.max(matrix)
    if max_val == min_val:
        return np.zeros(matrix.shape)
    return (matrix - min_val) / (max_val - min_val)

norm_gram_matrix = min_max_scale(gram_matrix)
norm_blosum_matrix = min_max_scale(blosum_matrix)

w_ctd = 0.4
w_blosum = 0.6
final_similarity_matrix = (w_ctd * norm_gram_matrix) + (w_blosum * norm_blosum_matrix)

distance_matrix = np.max(final_similarity_matrix) - final_similarity_matrix

output_filename = "final_distance_matrix.csv"
np.savetxt(output_filename, distance_matrix, delimiter=",")
print(f"Distance matrix successfully saved to {output_filename}")