import os
import sys
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from sklearn.decomposition import KernelPCA
from propy import Autocorrelation
import plotly.express as px
from sklearn.metrics.pairwise import rbf_kernel, pairwise_distances
from scipy.spatial import procrustes
import torch
from transformers import AutoTokenizer, AutoModel
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root) 

from app_components.config import COMPUTATION_DIR

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

def calculate_normalized_blosum_score(seq1: str, seq2: str) -> float:
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -1.0
    try:
        raw_score = aligner.score(seq1, seq2)
        min_len = min(len(seq1), len(seq2))
        return raw_score / min_len if min_len > 0 else 0.0
    except Exception:
        return 0.0

def min_max_scale(matrix):
    min_val, max_val = np.min(matrix), np.max(matrix)
    return (matrix - min_val) / (max_val - min_val) if max_val > min_val else np.zeros(matrix.shape)

def calculate_custom_kernel(df, gamma, lambda_w):
    print("-> Calculating Autocorrelation vectors...")
    df['autocorr_vector'] = df['h3_chain'].apply(
        lambda seq: compute_autocorr_vector(seq, charge550_dict, charge735_dict, hydro_dict)
    )
    df_filtered = df[df['autocorr_vector'].apply(len) > 0].reset_index(drop=True)
    feature_matrix = np.vstack(df_filtered['autocorr_vector'].values)

    print("-> Calculating RBF kernel matrix from autocorrelation vectors...")
    rbf_matrix = rbf_kernel(feature_matrix, gamma=gamma)

    print("-> Calculating pairwise normalized BLOSUM matrix...")
    num_sequences = len(df_filtered)
    blosum_matrix = np.zeros((num_sequences, num_sequences))
    sequences = df_filtered['h3_chain'].tolist()
    for i in range(num_sequences):
        for j in range(i, num_sequences):
            score = calculate_normalized_blosum_score(sequences[i], sequences[j])
            blosum_matrix[i, j] = blosum_matrix[j, i] = score

    print("-> Normalizing and combining matrices...")
    norm_rbf = min_max_scale(rbf_matrix)
    norm_blosum = min_max_scale(blosum_matrix)
    similarity_matrix = (lambda_w * norm_rbf) + ((1 - lambda_w) * norm_blosum)
    return similarity_matrix, df_filtered

def calculate_esm_kernel(df, model_name="facebook/esm2_t33_650M_UR50D"):
    print(f"-> Loading ESM model: {model_name}...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)

    print("-> Generating ESM embeddings...")
    sequences = df['h3_chain'].tolist()
    inputs = tokenizer(sequences, return_tensors="pt", padding=True, truncation=True, max_length=32)
    with torch.no_grad():
        outputs = model(**inputs)
    
    embeddings = outputs.last_hidden_state.mean(dim=1).numpy()
    
    print("-> Calculating pairwise cosine distances...")
    distance_matrix = pairwise_distances(embeddings, metric='cosine')
    
    print("-> Converting ESM distances to a similarity kernel...")
    similarity_matrix = 1 - distance_matrix
    return similarity_matrix, df

df_main = pd.read_csv(os.path.join(COMPUTATION_DIR, "cdr.csv"))
df_main.dropna(subset=['h3_chain'], inplace=True)
df_main = df_main.head(100)
df_main.reset_index(drop=True, inplace=True)

charge550_dict, charge735_dict, hydro_dict = prepare_charges_hydrophobicity()

gamma_val = 0.1
lambda_weight = 0.7

sim_custom, df_custom = calculate_custom_kernel(df_main.copy(), gamma=gamma_val, lambda_w=lambda_weight)
sim_esm, df_esm = calculate_esm_kernel(df_custom.copy())

print("-> Generating Kernel PCA projections for both methods...")
kpca = KernelPCA(n_components=2, kernel='precomputed', random_state=42)
embedding_custom_kpca = kpca.fit_transform(sim_custom)
embedding_esm_kpca = kpca.fit_transform(sim_esm)

print("-> Aligning Custom Kernel plot onto ESM benchmark plot using Procrustes...")
mtx1, mtx2, disparity = procrustes(embedding_esm_kpca, embedding_custom_kpca)
print(f"Disparity between point sets after alignment: {disparity:.4f}")

print("-> Creating linked scatter plot of aligned data...")
df_esm_aligned = pd.DataFrame({
    'x': mtx1[:, 0], 'y': mtx1[:, 1], 'method': 'ESM 650M (Reference)',
    'h3_chain': df_esm['h3_chain'], 'pdb_id': df_esm['pdb_id']
})
df_custom_aligned = pd.DataFrame({
    'x': mtx2[:, 0], 'y': mtx2[:, 1], 'method': 'Custom Kernel (Aligned)',
    'h3_chain': df_custom['h3_chain'], 'pdb_id': df_custom['pdb_id']
})
combined_plot_df = pd.concat([df_esm_aligned, df_custom_aligned])

fig = px.scatter(
    combined_plot_df, x='x', y='y', color='method',
    hover_data=['h3_chain', 'pdb_id'], title="Procrustes-Aligned Kernel PCA Projections: Custom vs ESM 650M"
)

for i in range(len(mtx1)):
    fig.add_shape(
        type='line',
        x0=mtx1[i, 0], y0=mtx1[i, 1],
        x1=mtx2[i, 0], y1=mtx2[i, 1],
        line=dict(color='rgba(128, 128, 128, 0.5)', width=1)
    )

fig.update_layout(xaxis_title="Principal Component 1", yaxis_title="Principal Component 2")
output_filename = "geary_procrustes_comparison.html"
fig.write_html(output_filename)
print(f"Aligned comparison plot saved successfully to '{output_filename}'.")