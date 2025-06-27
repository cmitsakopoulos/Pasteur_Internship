import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import umap
import hdbscan
import matplotlib.pyplot as plt
import os

st.set_page_config(layout="wide")
st.title("Advanced HDBSCAN Clustering with Parameter Estimation")

PATH_TO_DATA = "/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Internal_Files/h_chain_dist_matrix.npz"
EMBEDDING_FILE = 'saved_umap_embedding.csv'

@st.cache_data
def create_umap_embedding(path):
    try:
        with np.load(path) as data:
            if not data.files:
                return None, None
            matrix_key = data.files[0]
            distance_matrix = data[matrix_key]
        
        reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42)
        coords_2d = reducer.fit_transform(distance_matrix)
        
        df = pd.DataFrame({
            'x': coords_2d[:, 0],
            'y': coords_2d[:, 1],
            'id': [f'Point_{i}' for i in range(coords_2d.shape[0])]
        })
        return df, coords_2d
    except FileNotFoundError:
        st.error(f"Error: The file was not found at the path: {path}")
        return None, None
    except Exception as e:
        st.error(f"An error occurred: {e}")
        return None, None

if 'coords_df' not in st.session_state:
    st.session_state.coords_df = None
if 'coords_2d' not in st.session_state:
    st.session_state.coords_2d = None

st.sidebar.header("1. Load or Compute Data")
if os.path.exists(EMBEDDING_FILE):
    if st.sidebar.button("Load Saved UMAP Embedding"):
        with st.spinner("Loading from file..."):
            df = pd.read_csv(EMBEDDING_FILE)
            st.session_state.coords_df = df
            st.session_state.coords_2d = df[['x', 'y']].to_numpy()
            st.sidebar.success("Loaded saved embedding.")

if st.sidebar.button("Re-compute UMAP Embedding"):
    with st.spinner("Running UMAP... This may take a moment."):
        df, coords = create_umap_embedding(PATH_TO_DATA)
        if df is not None:
            st.session_state.coords_df = df
            st.session_state.coords_2d = coords
            df.to_csv(EMBEDDING_FILE, index=False)
            st.sidebar.success("UMAP complete and embedding saved.")

if st.session_state.coords_df is not None:
    coords_df = st.session_state.coords_df
    coords_2d = st.session_state.coords_2d

    st.sidebar.header("2. Parameter Estimation")
    with st.sidebar.expander("Explore Data Structure"):
        if st.button("Plot Condensed Cluster Tree"):
            clusterer = hdbscan.HDBSCAN(gen_min_span_tree=True).fit(coords_2d)
            fig, ax = plt.subplots(figsize=(10, 8))
            clusterer.condensed_tree_.plot(select_clusters=True, selection_palette=px.colors.qualitative.Plotly)
            st.pyplot(fig)

        if st.button("Calculate Cluster Stability (DBCV)"):
            if coords_2d is not None:
                sizes = range(5, 31, 2)
                scores = []
                with st.spinner("Calculating scores for different cluster sizes..."):
                    for size in sizes:
                        clusterer = hdbscan.HDBSCAN(min_cluster_size=size, gen_min_span_tree=True).fit(coords_2d)
                        score = clusterer.relative_validity_
                        scores.append(score)
            
                stability_df = pd.DataFrame({'min_cluster_size': sizes, 'DBCV_Score': scores})
                fig = px.line(stability_df, x='min_cluster_size', y='DBCV_Score', title='Cluster Stability (Higher is Better)', markers=True)
                st.plotly_chart(fig)
            else:
                st.warning("Data not loaded.")

    st.sidebar.header("3. Perform Final Clustering")
    min_cluster_size = st.sidebar.slider(
        'Minimum Cluster Size', min_value=2, max_value=50, value=10, step=1,
        help="The smallest size grouping that you wish to consider a cluster."
    )

    if st.sidebar.button("Run HDBSCAN Clustering"):
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
        st.session_state.cluster_labels = clusterer.fit_predict(coords_2d)
        labels = st.session_state.cluster_labels
        n_clusters_found = len(set(labels)) - (1 if -1 in labels else 0)
        noise_points = np.count_nonzero(labels == -1)
        st.success(f"Found {n_clusters_found} clusters and {noise_points} noise points.")

    if st.sidebar.button("Clear Clustering Results"):
        st.session_state.cluster_labels = None

    st.header("UMAP Projection with HDBSCAN Clustering")
    plot_df = coords_df.copy()
    color_map = None
    if 'cluster_labels' in st.session_state and st.session_state.cluster_labels is not None:
        plot_df['cluster'] = st.session_state.cluster_labels.astype(str)
        color_map = 'cluster'
    
    fig = px.scatter(
        plot_df, x='x', y='y', color=color_map,
        hover_data=['id', 'cluster'] if color_map else ['id'],
        color_discrete_map={"-1": "lightgrey"}
    )
    fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1))
    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Please load or compute a UMAP embedding to begin analysis.")