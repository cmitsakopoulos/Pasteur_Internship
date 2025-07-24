import streamlit as st
import pandas as pd
import plotly.express as px
from st_aggrid import AgGrid
import numpy as np
import hdbscan
from sklearn.cluster import SpectralClustering

from app_components.config import INTERNAL_FILES_DIR

def _perform_hdbscan(data_subset: pd.DataFrame, distance_matrix_subset: np.ndarray, min_cluster_size: int) -> pd.DataFrame:
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        metric='precomputed'
    )
    sub_cluster_labels = clusterer.fit_predict(distance_matrix_subset)
    
    result_df = data_subset.copy()
    result_df['sub_cluster_id'] = sub_cluster_labels
    return result_df

def _perform_spectral(data_subset: pd.DataFrame, affinity_matrix_subset: np.ndarray, n_clusters: int) -> pd.DataFrame:
    clusterer = SpectralClustering(
        n_clusters=n_clusters,
        affinity='precomputed',
        assign_labels='kmeans',
        random_state=69
    )
    sub_cluster_labels = clusterer.fit_predict(affinity_matrix_subset)
    
    result_df = data_subset.copy()
    result_df['sub_cluster_id'] = sub_cluster_labels
    return result_df

def render_sub_analysis_page():
    st.title("Protracted Cluster Analysis")

    method_map_files = {
        'Spectral Clustering': 'spectral_cluster',
        'HDBSCAN': 'snf_cluster'
    }
    selected_method_name = st.radio(
        "Select Overall Clustering Result:",
        options=list(method_map_files.keys()),
        horizontal=True,
    )
    
    results_basename = method_map_files[selected_method_name]
    html_path = INTERNAL_FILES_DIR  / f"{results_basename}_plot.html"
    try:
        with open(html_path, 'r', encoding='utf-8') as f:
            st.components.v1.html(f.read(), height=600, scrolling=True)
    except FileNotFoundError:
        st.warning(f"Static plot not found at: {html_path}")

    st.divider()
    st.header("Recursive Sub-Cluster Analysis")
    st.markdown("Select a parent cluster, then choose a method to re-cluster its members.")
    
    parquet_path = INTERNAL_FILES_DIR  / f"{results_basename}.parquet"
    matrix_path = INTERNAL_FILES_DIR / "fused_affinity_matrix.npz"

    try:
        clustered_df = pd.read_parquet(parquet_path)
        affinity_matrix_full = np.load(matrix_path)['matrix']
    except FileNotFoundError:
        st.error(f"Missing required file: {parquet_path} or {matrix_path}")
        st.stop()

    clustered_df['cluster_id'] = clustered_df['cluster_id'].astype(str)
    
    cluster_ids = sorted([c for c in clustered_df['cluster_id'].unique()])
    display_map = {c: f"Cluster {c}" for c in cluster_ids if c != '-1'}
    display_map['-1'] = "Noise (-1)"
    display_options = sorted(display_map.values(), key=lambda x: int(x.split(' ')[-1].strip('()')) if 'Cluster' in x else -1)

    if not cluster_ids:
        st.warning("No analysable clusters found in the data file.")
        st.stop()
    
    col1, col2 = st.columns(2)
    with col1:
        selected_display_name = st.selectbox("Select a Parent Group to Analyse:", display_options)
        
    reverse_display_map = {v: k for k, v in display_map.items()}
    selected_id = reverse_display_map[selected_display_name]

    if selected_id == '-1':
        st.info("Noise points represent outliers and cannot be sub-clustered.")
    else:
        with col2:
            sub_method = st.radio("Sub-Clustering Method:", ["HDBSCAN", "Spectral Clustering"], horizontal=True)

        if sub_method == "HDBSCAN":
            min_size = st.slider("HDBSCAN: Minimum Cluster Size", min_value=2, max_value=50, value=5)
        else:
            k_clusters = st.slider("Spectral: Number of Sub-Clusters (k)", min_value=2, max_value=10, value=3)
        
        if st.button(f"Analyse Sub-Cluster {selected_id}", type="primary", use_container_width=True):
            st.session_state.sub_analysis_params = {
                "cluster_id": selected_id,
                "method": sub_method,
                "params": {"min_cluster_size": min_size} if sub_method == "HDBSCAN" else {"n_clusters": k_clusters}
            }
            if 'sub_analysis_result' in st.session_state:
                del st.session_state['sub_analysis_result']

    if 'sub_analysis_params' in st.session_state and st.session_state.sub_analysis_params['cluster_id'] != '-1':
        if 'sub_analysis_result' not in st.session_state:
            params = st.session_state.sub_analysis_params
            parent_cluster_id = params['cluster_id']
            method = params['method']

            with st.spinner(f"Running {method} on members of Cluster {parent_cluster_id}..."):
                data_subset = clustered_df[clustered_df['cluster_id'] == parent_cluster_id].copy()
                subset_indices = data_subset.index.to_numpy()
                affinity_matrix_subset = affinity_matrix_full[subset_indices, :][:, subset_indices]

                if method == "HDBSCAN":
                    distance_matrix_subset = 1 - affinity_matrix_subset
                    result_df = _perform_hdbscan(data_subset, distance_matrix_subset, **params['params'])
                else:
                    result_df = _perform_spectral(data_subset, affinity_matrix_subset, **params['params'])
                
                st.session_state.sub_analysis_result = result_df
                st.session_state.sub_analysis_display_title = f"Sub-Analysis of Cluster {parent_cluster_id} using {method}"

    if 'sub_analysis_result' in st.session_state:
        st.subheader(st.session_state.sub_analysis_display_title)
        result_df = st.session_state.sub_analysis_result
        
        result_df['sub_cluster_id'] = result_df['sub_cluster_id'].astype(str).replace({'-1': 'Noise'})

        fig = px.scatter(
            result_df, x='umap_x', y='umap_y', color='sub_cluster_id',
            color_discrete_map={'Noise': 'lightgrey'},
            hover_name='h3_chain', title=f"Interactive Plot of New Sub-Clusters"
        )
        fig.update_traces(marker=dict(size=8, line=dict(width=1, color='Black')))
        fig.update_layout(legend_title_text='Sub-Cluster ID')
        st.plotly_chart(fig, use_container_width=True)

        st.write("Sub-Clustering Results:")
        all_columns = result_df.columns.tolist()
        default_cols = [col for col in ['sub_cluster_id', 'id', 'pdb_id', 'h3_chain'] if col in all_columns]
        
        selected_columns = st.multiselect(
            "Select metadata columns to display:",
            options=all_columns,
            default=default_cols
        )
        if selected_columns:
            AgGrid(result_df[selected_columns].fillna("NULL"), height=300, fit_columns_on_grid_load=True)
        else:
            st.warning("Please select at least one column to display.")