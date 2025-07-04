import numpy as np
import pandas as pd
import umap
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from propy import Autocorrelation
import plotly.express as px

def prepare_charges_hydrophobicity():
    aas = "ACDEFGHIKLMNPQRSTVWY"
    aa_list = list(aas)
    net_charge_at_pH735 = {aa: ProteinAnalysis(aa).charge_at_pH(7.35) for aa in aa_list}
    net_charge_at_pH550 = {aa: ProteinAnalysis(aa).charge_at_pH(5.5) for aa in aa_list}
    hydrophobicity_scale = ProtParamData.kd
    return net_charge_at_pH550, net_charge_at_pH735, hydrophobicity_scale

def compute_autocorr_vector(sequence: str, charge55, charge735, hydro) -> np.ndarray:
    try:
        if not isinstance(sequence, str) or len(sequence) < 2:
            return np.array([])
        
        ph550_vals = list(Autocorrelation.CalculateEachGearyAuto(sequence, charge55, "net_charges").values())
        ph735_vals = list(Autocorrelation.CalculateEachGearyAuto(sequence, charge735, "net_charges").values())
        hydro_vals = list(Autocorrelation.CalculateEachGearyAuto(sequence, hydro, "hydrophobicity").values())
        
        return np.concatenate([ph550_vals, ph735_vals, hydro_vals])
    except Exception:
        return np.array([])

def calculate_blosum_score(seq1: str, seq2: str) -> float:
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -1.0
    try:
        return aligner.score(seq1, seq2)
    except Exception:
        return 0.0

def min_max_scale(matrix):
    min_val, max_val = np.min(matrix), np.max(matrix)
    return (matrix - min_val) / (max_val - min_val) if max_val > min_val else np.zeros(matrix.shape)

df_main = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/cdr.csv")

charge550_dict, charge735_dict, hydro_dict = prepare_charges_hydrophobicity()

print("-> Calculating Autocorrelation vectors...")
df_main['autocorr_vector'] = df_main['h3_chain'].apply(
    lambda seq: compute_autocorr_vector(seq, charge550_dict, charge735_dict, hydro_dict)
)
df_filtered = df_main[df_main['autocorr_vector'].apply(len) > 0].reset_index(drop=True)

feature_matrix = np.vstack(df_filtered['autocorr_vector'].values)
print("-> Calculating Gram matrix (dot products)...")
gram_matrix = feature_matrix @ feature_matrix.T

print("-> Calculating pairwise BLOSUM matrix...")
num_sequences = len(df_filtered)
blosum_matrix = np.zeros((num_sequences, num_sequences))
sequences = df_filtered['h3_chain'].tolist()
for i in range(num_sequences):
    for j in range(i, num_sequences):
        score = calculate_blosum_score(sequences[i], sequences[j])
        blosum_matrix[i, j] = blosum_matrix[j, i] = score

print("-> Normalizing and combining matrices...")
norm_gram = min_max_scale(gram_matrix)
norm_blosum = min_max_scale(blosum_matrix)
similarity_matrix = (0.5 * norm_gram) + (0.5 * norm_blosum)

print("-> Converting similarity matrix to distance matrix for UMAP...")
distance_matrix = np.max(similarity_matrix) - similarity_matrix

print("-> Generating projection using UMAP...")
reducer = umap.UMAP(
    n_components=2,
    metric='precomputed',
    random_state=42
)
embedding_umap = reducer.fit_transform(distance_matrix)

plot_df = pd.DataFrame({
    'x': embedding_umap[:, 0], 'y': embedding_umap[:, 1],
    'pdb_id': df_filtered['pdb_id'], 'h3_chain': df_filtered['h3_chain']
})
fig = px.scatter(
    plot_df, x='x', y='y', hover_data=['pdb_id', 'h3_chain'],
    title="UMAP Projection of Autocorrelation + BLOSUM Similarity"
)

output_filename = "autocorrelation_kernel_umap_plot.html"
fig.write_html(output_filename)
print(f"Plot saved successfully to '{output_filename}'.")