import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from sklearn.manifold import MDS
from sklearn.cluster import AgglomerativeClustering

st.set_page_config(layout="wide")
st.title("Interactive Hierarchical Clustering of Sequence Data")

PATH_TO_DATA = "/Users/chrismitsacopoulos/Desktop/Pasteur_Internship/Internal_Files/h3_chain_dist_matrix.npz"

@st.cache_data
def create_mds_coordinates(path):
    try:
        with np.load(path) as data:
            if not data.files:
                return None
            matrix_key = data.files[0]
            distance_matrix = data[matrix_key]

        mds = MDS(
            n_components=2,
            metric=False,
            n_init=10,
            dissimilarity='precomputed',
            random_state=42
        )
        coords_2d = mds.fit_transform(distance_matrix)

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

coords_df = create_mds_coordinates(PATH_TO_DATA)

if 'cluster_labels' not in st.session_state:
    st.session_state.cluster_labels = None

st.sidebar.header("Clustering Controls")
n_clusters = st.sidebar.number_input(
    'Select number of clusters to extract',
    min_value=2,
    max_value=20,
    value=5,
    step=1
)

if st.sidebar.button("Run Agglomerative Clustering"):
    if coords_df is not None:
        agg_clustering = AgglomerativeClustering(n_clusters=n_clusters)
        st.session_state.cluster_labels = agg_clustering.fit_predict(coords_df[['x', 'y']])
        st.success(f"Successfully clustered data into {n_clusters} groups.")
    else:
        st.sidebar.error("Data not loaded. Cannot run clustering.")

if st.sidebar.button("Clear Clustering Results"):
    st.session_state.cluster_labels = None

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
        title="MDS Plot with Agglomerative Clustering"
    )

    fig.update_traces(marker=dict(size=8))
    fig.update_layout(
        yaxis=dict(scaleanchor="x", scaleratio=1),
        xaxis_title="MDS Dimension 1",
        yaxis_title="MDS Dimension 2"
    )

    st.plotly_chart(fig, use_container_width=True)