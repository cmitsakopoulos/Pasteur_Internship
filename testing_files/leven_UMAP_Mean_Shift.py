import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import umap
from sklearn.cluster import MeanShift

st.set_page_config(layout="wide")
st.title("Interactive Clustering on a UMAP Projection")

PATH_TO_DATA = "/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Internal_Files/h3_chain_dist_matrix.npz"

@st.cache_data
def create_umap_embedding(path):
    try:
        with np.load(path) as data:
            if not data.files:
                return None
            matrix_key = data.files[0]
            distance_matrix = data[matrix_key]

        st.write("Creating 2D embedding using UMAP...")
        reducer = umap.UMAP(
            n_components=2,
            metric='precomputed',
            random_state=42
        )
        coords_2d = reducer.fit_transform(distance_matrix)
        st.write("UMAP embedding complete.")

        num_points = distance_matrix.shape[0]
        df = pd.DataFrame({
            'x': coords_2d[:, 0],
            'y': coords_2d[:, 1],
            'id': [f'Point_{i}' for i in range(num_points)]
        })
        return df
    except FileNotFoundError:
        st.error(f"Error: The file was not found at the path: {path}")
        return None
    except Exception as e:
        st.error(f"An error occurred: {e}")
        return None

coords_df = create_umap_embedding(PATH_TO_DATA)

if 'cluster_labels' not in st.session_state:
    st.session_state.cluster_labels = None
    st.session_state.n_clusters_found = 0

st.sidebar.header("Clustering Controls")
st.sidebar.write("Mean Shift automatically discovers the number of clusters.")

if st.sidebar.button("Run Mean Shift Clustering"):
    if coords_df is not None:
        mean_shift = MeanShift(bin_seeding=True)
        st.session_state.cluster_labels = mean_shift.fit_predict(coords_df[['x', 'y']])
        n_clusters_found = len(np.unique(st.session_state.cluster_labels))
        st.session_state.n_clusters_found = n_clusters_found
        st.success(f"Clustering found {n_clusters_found} groups.")
    else:
        st.sidebar.error("Data not loaded. Cannot run clustering.")

if st.sidebar.button("Clear Clustering Results"):
    st.session_state.cluster_labels = None
    st.session_state.n_clusters_found = 0

if coords_df is not None:
    plot_df = coords_df.copy()
    color_map = None

    if st.session_state.cluster_labels is not None:
        plot_df['cluster'] = st.session_state.cluster_labels.astype(str)
        color_map = 'cluster'

    fig = px.scatter(
        plot_df,
        x='x',
        y='y',
        color=color_map,
        hover_data=['id', 'cluster'] if color_map else ['id'],
        title="UMAP Projection with Mean Shift Clustering"
    )

    fig.update_traces(marker=dict(size=8))
    fig.update_layout(
        yaxis=dict(scaleanchor="x", scaleratio=1),
        xaxis_title="UMAP Dimension 1",
        yaxis_title="UMAP Dimension 2"
    )

    st.plotly_chart(fig, use_container_width=True)