import numpy as np
import pandas as pd
import umap
import plotly.express as px
from scipy.spatial import procrustes
import torch

from propy.PyPro import GetProDes
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from sklearn.metrics.pairwise import pairwise_distances
from transformers import AutoTokenizer, AutoModel

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
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    try:
        return aligner.score(seq1, seq2)
    except Exception:
        return 0.0

def min_max_scale(matrix):
    min_val, max_val = np.min(matrix), np.max(matrix)
    return (matrix - min_val) / (max_val - min_val) if max_val > min_val else np.zeros(matrix.shape)

def calculate_custom_kernel_distance(df):
    print("-> Calculating CTD vectors...")
    df['ctd_vector'] = df['h3_chain'].apply(compute_ctd_vector)
    df_filtered = df[df['ctd_vector'].apply(len) > 0].reset_index(drop=True)

    feature_matrix = np.vstack(df_filtered['ctd_vector'].values)
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

    similarity_matrix = (1 * norm_gram) + (1 * norm_blosum)
    distance_matrix = np.max(similarity_matrix) - similarity_matrix
    return distance_matrix, df_filtered

def calculate_esm_distance(df, model_name="facebook/esm2_t12_35M_UR50D"):
    print(f"-> Loading ESM model: {model_name}...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)

    print("-> Generating ESM embeddings...")
    sequences = df['h3_chain'].tolist()
    inputs = tokenizer(sequences, return_tensors="pt", padding=True, truncation=True)
    with torch.no_grad():
        outputs = model(**inputs)

    embeddings = outputs.last_hidden_state.mean(dim=1).numpy()

    print("-> Calculating pairwise cosine distances...")
    distance_matrix = pairwise_distances(embeddings, metric='cosine')
    return distance_matrix, df

df_main = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/cdr.csv")
df_main.dropna(subset=['h3_chain'], inplace=True)
df_main = df_main.head(100)
df_main.reset_index(drop=True, inplace=True)

dist_custom, df_custom = calculate_custom_kernel_distance(df_main.copy())
dist_esm, df_esm = calculate_esm_distance(df_main.copy())

print("-> Generating UMAP projections...")
reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42)
embedding_custom = reducer.fit_transform(dist_custom)
embedding_esm = reducer.fit_transform(dist_esm)

print("-> Aligning plots using Procrustes Analysis...")
mtx1, mtx2, disparity = procrustes(embedding_custom, embedding_esm)
print(f"Disparity between point sets after alignment: {disparity:.4f}")

print("-> Creating linked scatter plot of aligned data...")
df_custom_aligned = pd.DataFrame({
    'x': mtx1[:, 0], 'y': mtx1[:, 1], 'method': 'Custom Kernel (Reference)',
    'h3_chain': df_custom['h3_chain'], 'pdb_id': df_custom['pdb_id']
})
df_esm_aligned = pd.DataFrame({
    'x': mtx2[:, 0], 'y': mtx2[:, 1], 'method': 'ESM (Aligned)',
    'h3_chain': df_esm['h3_chain'], 'pdb_id': df_esm['pdb_id']
})
combined_plot_df = pd.concat([df_custom_aligned, df_esm_aligned])

fig = px.scatter(
    combined_plot_df, x='x', y='y', color='method',
    hover_data=['h3_chain', 'pdb_id'], title="Procrustes-Aligned UMAP Projections"
)

for i in range(len(mtx1)):
    fig.add_shape(
        type='line',
        x0=mtx1[i, 0], y0=mtx1[i, 1],
        x1=mtx2[i, 0], y1=mtx2[i, 1],
        line=dict(color='rgba(128, 128, 128, 0.5)', width=1)
    )

fig.update_layout(xaxis_title="Dimension 1", yaxis_title="Dimension 2")
output_filename = "procrustes_comparison_plot_650M.html"
fig.write_html(output_filename)
print(f"Aligned comparison plot saved successfully to '{output_filename}'.")