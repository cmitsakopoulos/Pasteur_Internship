import os
import sys
import numpy as np
import pandas as pd
import ast
from Bio.Align import PairwiseAligner, substitution_matrices
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler
import snf
import umap
import plotly.express as px

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

from app_components.config import COMPUTATION_DIR

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
df_main.dropna(subset=['h3_chain', 'h3_ctd'], inplace=True)
df_main = df_main.head(100).reset_index(drop=True)

print("-> Parsing pre-computed CTD vectors...")
df_main['ctd_vector'] = df_main['h3_ctd'].apply(ast.literal_eval)
ctd_feature_matrix = np.array(df_main['ctd_vector'].tolist())

print("-> Standardizing CTD features...")
scaler = StandardScaler()
scaled_ctd_features = scaler.fit_transform(ctd_feature_matrix)

print("-> Computing distance matrices...")
ctd_distance_matrix = pairwise_distances(scaled_ctd_features, metric='manhattan')

sequences = df_main['h3_chain'].tolist()
num_sequences = len(sequences)
blosum_similarity_matrix = np.zeros((num_sequences, num_sequences))
for i in range(num_sequences):
    for j in range(i, num_sequences):
        score = calculate_blosum_score(sequences[i], sequences[j])
        blosum_similarity_matrix[i, j] = blosum_similarity_matrix[j, i] = score

blosum_distance_matrix = np.max(blosum_similarity_matrix) - blosum_similarity_matrix

print("-> Creating affinity matrices from distance matrices...")
K_NEIGHBOURS = 20
MU_PARAM = 0.5

affinity_ctd = snf.compute.affinity_matrix(ctd_distance_matrix, K=K_NEIGHBOURS, mu=MU_PARAM)
affinity_blosum = snf.compute.affinity_matrix(blosum_distance_matrix, K=K_NEIGHBOURS, mu=MU_PARAM)

affinity_matrices = [affinity_blosum, affinity_ctd]

print("-> Fusing networks using SNF...")
fused_affinity_matrix = snf.snf(affinity_matrices, K=K_NEIGHBOURS)

print("-> Projecting fused affinity matrix into 2D space using UMAP...")
fused_distance_matrix = 1 - fused_affinity_matrix
umap_reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42, n_neighbors=15, min_dist=0.1)
embedding_2d = umap_reducer.fit_transform(fused_distance_matrix)

print("-> Generating interactive plot...")
plot_df = pd.DataFrame(embedding_2d, columns=['UMAP1', 'UMAP2'])
plot_df['h3_chain'] = sequences
plot_df['pdb_id'] = df_main['pdb_id']

fig = px.scatter(
    plot_df,
    x='UMAP1',
    y='UMAP2',
    hover_data=['h3_chain', 'pdb_id'],
    title='UMAP Projection of SNF Fused Similarity Network'
)
fig.update_traces(marker=dict(size=8))
fig.update_layout(
    xaxis_title="UMAP Component 1",
    yaxis_title="UMAP Component 2"
)

output_filename = "snf_projection_fixed.html"
fig.write_html(output_filename)
print(f"\n--- Computation Complete ---")
print(f"Interactive plot saved successfully to '{output_filename}'.")