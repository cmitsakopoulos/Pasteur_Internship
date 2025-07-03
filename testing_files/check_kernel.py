import umap
import numpy as np
import plotly.express as px
import pandas as pd

df = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Computation_Deposit/cdr.csv")

distance_df = pd.read_csv("/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/final_distance_matrix.csv", header=None)
distance_matrix = distance_df.to_numpy()
print("Generating UMAP projection for visualisation...")
reducer = umap.UMAP(
    n_components=2,
    metric='precomputed',
    random_state=42
)

embedding = reducer.fit_transform(distance_matrix)

plot_df = pd.DataFrame({
    'x': embedding[:, 0],
    'y': embedding[:, 1],
    'h3_chain': df["h3_chain"]
})

fig = px.scatter(
    plot_df,
    x='x',
    y='y',
    hover_data=['h3_chain'],
    title="UMAP Projection of Custom Kernel Distances"
)

output_filename = "umap_projection.html"
fig.write_html(output_filename)

print(f"Plot saved successfully to '{output_filename}'. You can open this file in your browser.")