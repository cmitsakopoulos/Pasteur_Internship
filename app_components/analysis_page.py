import streamlit as st
import pandas as pd
import time
import os
import threading
from io import StringIO
import contextlib
from streamlit_autorefresh import st_autorefresh

from app_components.config import IMAGES_DIR, INTERNAL_FILES_DIR

from app_components.database_pipeline import (
    FileLogger, Pipeline, Concatenation, CDRComputation, 
    AntigenComputation, FlattenDuplicates, Write, RmPurificationTags, 
    AssignIDs, ComputeRelationships, CleanUp, 
    PreWalker, Walker, HDBSCAN, ComputeDistanceMatrices, TuneSNFParameters, FuseAndProject, Spectral_Clustering, MantelTest, Procrustes_Analysis
)
from app_components.function_dump import inspect_summary, inspect_verbose

RECIPES = {
    "Standard": [CleanUp, PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation, AssignIDs, ComputeRelationships, Write],
    "Rerun": [PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
    "Distance Matrix": [ComputeDistanceMatrices],
    "Generate Scatter Plot": [ComputeDistanceMatrices, TuneSNFParameters, FuseAndProject],
    "Mantel Test": [MantelTest],
    "Procrustes": [Procrustes_Analysis],
    "HDBSCAN": [HDBSCAN],
    "SpectralClustering": [Spectral_Clustering],
}

def build_pipeline(recipe: str, path: str, logger: FileLogger, params: dict = None) -> Pipeline:
    if params is None:
        params = {}
    logger.log(f"Building '{recipe}' pipeline for path: {path} with params: {params}")
    steps = []
    for step_class in RECIPES[recipe]:
        if step_class is PreWalker:
            steps.append(step_class(input_path=path, logger=logger))
        elif step_class is HDBSCAN:
            min_size = params.get('min_cluster_size', 15)
            steps.append(step_class(logger=logger, min_cluster_size=min_size))
        elif step_class is Spectral_Clustering:
            desired_n = params.get('n_clusters', 5)
            steps.append(step_class(logger=logger, n_clusters=desired_n))
        else:
            steps.append(step_class(logger=logger))
    return Pipeline(steps, logger)

def run_pipeline_worker(job_id, recipe, path, params: dict = None):
    logger = None
    try:
        logger = FileLogger(job_id)
        pipeline = build_pipeline(recipe, path, logger, params)
        pipeline.run()
    except BaseException as e:
        if logger:
            logger.log(f"\n--- A FATAL ERROR OCCURRED IN WORKER ---\n{type(e).__name__}\n{str(e)}")
    finally:
        if logger:
            logger.close()

def get_subdirectories(path='.'):
    try:
        dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
        return sorted(dirs)
    except FileNotFoundError:
        return []

def render_dashboard_page():

    st_autorefresh(interval=2000, key="global_refresher")
    if "text_output" not in st.session_state: st.session_state.text_output = ""
    if "running_job_id" not in st.session_state: st.session_state.running_job_id = None
    if "active_thread" not in st.session_state: st.session_state.active_thread = None

    control_col, dashboard_col = st.columns((1, 2))

    with control_col:
        st.image(str(IMAGES_DIR / "logo.png"), width=230)
        st.title("Control Panel")

        if st.session_state.running_job_id and st.session_state.active_thread:
            if not st.session_state.active_thread.is_alive():
                job_id = st.session_state.running_job_id
                log_filename = INTERNAL_FILES_DIR / f"{job_id}.log"

                if os.path.exists(log_filename):
                    try:
                        with open(log_filename, 'r', encoding='utf-8') as f:
                            st.session_state.text_output = f.read()
                        os.remove(log_filename)
                        st.toast("Pipeline run finished!", icon="🕵️‍♂️")
                    except Exception as e:
                        st.error(f"Could not read or delete log file: {e}")
                else:
                    st.session_state.text_output = f"Job finished, but log file ('{log_filename}') was not found."

                st.session_state.running_job_id = None
                st.session_state.active_thread = None
        
        is_running = bool(st.session_state.running_job_id) 

        with st.expander("**1. Process Raw Data**", expanded=False):
            recipe = st.radio("Select Recipe", ('Standard', 'Rerun'),
                captions=["Clears prior results and runs fresh.", "Runs on new files only."],
                horizontal=True)
            path_run = st.text_input("Input Directory", value=str(INTERNAL_FILES_DIR),
                help="Specify the directory containing the raw data to be processed.")

            if st.button("Execute Pipeline", type="primary", use_container_width=True, disabled=is_running):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id
                    clean_path = path_run.strip().strip('"') #learned this the hard way
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, recipe, clean_path))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info(f"Started '{recipe}' pipeline in the background.")
                else:
                    st.warning("Another pipeline is already running.")

        with st.expander("**2. Separate CTD/BLOSUM Descriptor Analysis**", expanded=False):
            path_analysis = str(INTERNAL_FILES_DIR)

            st.subheader("Distance Matrices for Downstream Analysis")
            if st.button("Calculate Distance Matrices", key= "dist_separate", use_container_width=True, disabled=is_running):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Distance Matrix", path_analysis))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info("Started distance matrix calculation.")
                else: st.warning("Another pipeline is already running.")

            st.divider()
            st.subheader("Mantel Correlation Test Between Distance Matrices")
            if st.button("Perform Mantel Comparison", use_container_width=True, disabled=is_running):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Mantel Test", path_analysis))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info("Started Mantel test.")
                else: st.warning("Another pipeline is already running.")

            st.divider()
            st.subheader("Visual-Procrustes-Comparison on MDS Projection of Distance Matrices")
            if st.button("Perform Procrustes Comparison", use_container_width=True, disabled=is_running):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Procrustes", path_analysis))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info("Started Procrustes analysis.")
                else: st.warning("Another pipeline is already running.")


        with st.expander("**3. Similarity Network Fusion CTD+BLOSUM Analysis**", expanded=False):
            path_analysis = str(INTERNAL_FILES_DIR)

            st.subheader("Distance and Approximate Space Projection")
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Calculate Distance Matrices", key= "dist_SNF", use_container_width=True, disabled=is_running):
                    if not st.session_state.running_job_id:
                        st.session_state.text_output = ""
                        job_id = f"job_{time.time()}"
                        st.session_state.running_job_id = job_id
                        thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Distance Matrix", path_analysis))
                        st.session_state.active_thread = thread
                        thread.start()
                        st.info("Started distance matrix calculation.")
                    else: st.warning("Another pipeline is already running.")
            with col2:
                if st.button("Generate UMAP Plot", use_container_width=True, disabled=is_running):
                    if not st.session_state.running_job_id:
                        st.session_state.text_output = ""
                        job_id = f"job_{time.time()}"
                        st.session_state.running_job_id = job_id
                        thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Generate Scatter Plot", path_analysis))
                        st.session_state.active_thread = thread
                        thread.start()
                        st.info("Started UMAP generation.")
                    else: st.warning("Another pipeline is already running.")

            st.divider()
            st.subheader("HDBSCAN Clustering")
            min_size = st.number_input("Minimum Cluster Size", min_value=2, value=15, step=1, help="Set minimum cluster size for HDBSCAN.")
            
            if st.button("Perform HDBSCAN Clustering", type="primary", use_container_width=True, disabled=is_running):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id 
                    params = {'min_cluster_size': min_size}
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "HDBSCAN", path_analysis, params))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info("Started HDBSCAN clustering.")
                else: st.warning("Another pipeline is already running.")

            st.divider()
            st.subheader("Spectral Clustering")
            desired_n = st.number_input("Desired Clusters", min_value=2, value=7, step=1, help="Set amount of desired clusters.")
            
            if st.button("Perform Spectral Clustering", type="primary", use_container_width=True, disabled=is_running):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id 
                    params = {'n_clusters': desired_n}
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "SpectralClustering", path_analysis, params))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info("Started Spectral clustering.")
                else: st.warning("Another pipeline is already running.")

        with st.expander("**4. Inspect a Data File**"):
            uploaded_file = st.file_uploader("Upload a data file for analysis", type=["csv"], help="Upload a processed or external CSV file.")
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Quick Summary", use_container_width=True, disabled=not uploaded_file):
                    st.session_state.text_output = ""
                    log_capture_string = StringIO()
                    try:
                        with contextlib.redirect_stdout(log_capture_string):
                            df = pd.read_csv(uploaded_file)
                            inspect_summary(df)
                        st.session_state.text_output = log_capture_string.getvalue()
                        st.toast("Summary in Logs", icon="🕵️‍♂️")
                    except Exception as e: st.error(f"Failed to generate summary: {e}")
            with col2:
                if st.button("Detailed Analysis", use_container_width=True, disabled=not uploaded_file):
                    st.session_state.text_output = ""
                    log_capture_string = StringIO()
                    try:
                        with contextlib.redirect_stdout(log_capture_string):
                            df = pd.read_csv(uploaded_file)
                            inspect_verbose(df)
                        st.session_state.text_output = log_capture_string.getvalue()
                        st.toast("Expanded View in Logs", icon="🕵️‍♂️")
                    except Exception as e: st.error(f"Failed to run analysis: {e}")
            if not uploaded_file: st.caption("Upload a file to enable inspection options.")

    with dashboard_col:

        st.header("Activity & Results")

        if st.session_state.running_job_id:
                spinner_css = """
                    <style>
                    .spinner {
                        border: 4px solid rgba(255, 255, 255, 0.3);
                        border-left-color: #09ab3b; 
                        border-radius: 50%;
                        width: 30px;
                        height: 30px;
                        animation: spin 1s linear infinite;
                    }
                    @keyframes spin {
                        to {
                            transform: rotate(360deg);
                        }
                    }
                    </style>
                    """
                spinner_html = "<div class='spinner'></div>"

                spinner_col, text_col = st.columns([1, 10])
                with spinner_col:
                    st.markdown(spinner_css + spinner_html, unsafe_allow_html=True)
                with text_col:
                    st.info("Pipeline is running in the background. This page will auto-refresh.")
        

        tab_plots, tab_logs = st.tabs(["Visualisations", "Logs"])
        
        with tab_plots:

            PLOTS = {
                "UMAP Projection": {"file": "snf_umap_projection.html", "header": "UMAP (Non-Linear) Projection of Relative Distance"},
                "HDBSCAN Clustering": {"file": "snf_cluster_plot.html", "header": "UMAP (Non-Linear) Projection on HDBSCAN Clustering of Distance Matrices"},
                "Spectral Clustering": {"file": "spectral_cluster_plot.html", "header": "UMAP (Non-Linear) Projection on Spectral Clustering of Distance Matrices"},
                "Procrustes Comparison": {"file": "procrustes_comparison.html", "header": "MDS (Linear) Projection and Procrustes"},
                "Mantel Test": {"file": "mantel_test_plot.html", "header": "Mantel Stastical Correlation Between Distance Matrices"}
            }

            selected_plot_name = st.selectbox(
                "Choose a visualisation to display:",
                options=list(PLOTS.keys())
            )

            if selected_plot_name:
                plot_info = PLOTS[selected_plot_name]
                plot_file_path = INTERNAL_FILES_DIR / plot_info["file"]
                
                if os.path.exists(plot_file_path):
                    st.subheader(plot_info["header"])
                    with open(plot_file_path, 'r', encoding='utf-8') as f:
                        html_content = f.read()
                    st.components.v1.html(html_content, height=550, scrolling=True)
                else:
                    st.warning(f"Plot '{plot_info['header']}' has not been generated yet. Please run the corresponding step in the 'Analysis & Visualisation' panel. Each graph in the drop down has its own step with which they are produced for stability purposes.")
        with tab_logs:
            st.subheader("Operations Log")
            if st.session_state.text_output:
                st.code(st.session_state.text_output, language='bash')
            elif not st.session_state.running_job_id:
                st.info("Output from a pipeline run or file inspection will be shown here.")
            else:
                st.info("Job is in progress...")