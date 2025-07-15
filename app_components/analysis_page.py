import streamlit as st
import pandas as pd
import time
import psutil
import os
import threading
from io import StringIO
import contextlib
from streamlit_autorefresh import st_autorefresh
from streamlit_extras.metric_cards import style_metric_cards
from streamlit_extras.colored_header import colored_header
from app_components.config import IMAGES_DIR, INTERNAL_FILES_DIR

from app_components.database_pipeline import (
    FileLogger, Pipeline, Concatenation, CDRComputation, 
    AntigenComputation, FlattenDuplicates, Write, RmPurificationTags, 
    AssignIDs, ComputeRelationships, CleanUp, 
    PreWalker, Walker, GetDensity, Dendrogram, HDBSCAN, ComputeDistanceMatrices, TuneSNFParameters, FuseAndProject
)
from app_components.function_dump import inspect_summary, inspect_verbose

RECIPES = {
    "Standard": [CleanUp, PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation, AssignIDs, ComputeRelationships, Write],
    "Rerun": [PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
    "Distance Matrix": [ComputeDistanceMatrices],
    "Generate Scatter Plot": [ComputeDistanceMatrices, TuneSNFParameters, FuseAndProject],
    "Identify MinClusterNum": [GetDensity],
    "Dendrogram": [Dendrogram],
    "Perform Clustering": [HDBSCAN]
}

def build_pipeline(recipe: str, path: str, logger: FileLogger, params: dict = None) -> Pipeline:
    if params is None:
        params = {}
    logger.log(f"Building '{recipe}' pipeline for path: {path} with params: {params}")
    steps = []
    for step_class in RECIPES[recipe]:
        if step_class is PreWalker:
            steps.append(step_class(input_path=path, logger=logger))
        elif step_class is Dendrogram or step_class is HDBSCAN:
            min_size = params.get('min_cluster_size', 15)
            steps.append(step_class(logger=logger, min_cluster_size=min_size))
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

    if "text_output" not in st.session_state: st.session_state.text_output = ""
    if "running_job_id" not in st.session_state: st.session_state.running_job_id = None
    if "active_thread" not in st.session_state: st.session_state.active_thread = None

    control_col, dashboard_col = st.columns((1, 2))

    with control_col:
        st.image("./Images/logo.png", width=230)
        st.title("Control Panel")

        if st.session_state.running_job_id and st.session_state.active_thread:
            if not st.session_state.active_thread.is_alive():
                job_id = st.session_state.running_job_id
                log_filename = os.path.abspath(os.path.join("Internal_Files", f"{job_id}.log"))

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
        
        with st.expander("**2. Analysis**", expanded=False):
            path_analysis = str(INTERNAL_FILES_DIR)

            st.subheader("Relative Distance and Space Projection")
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Calculate Distance Matrix", use_container_width=True, disabled=is_running):
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
            st.subheader("Density Analysis-Hyperparametre Tuning")
            
            if st.button("Analyse Cluster Stability", use_container_width=True, disabled=is_running):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Identify MinClusterNum", path_analysis))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info("Started stability analysis.")
                else: st.warning("Another pipeline is already running.")

            st.divider()
            st.subheader("HDBSCAN Clustering")
            min_size = st.number_input("Minimum Cluster Size", min_value=2, value=15, step=1, help="Set the minimum cluster size for Dendrogram and HDBSCAN steps.")
            
            col3 = st.columns(1)
            with col3:
                if st.button("Perform HDBSCAN Clustering", type="primary", use_container_width=True, disabled=is_running):
                    if not st.session_state.running_job_id:
                        st.session_state.text_output = ""
                        job_id = f"job_{time.time()}"
                        st.session_state.running_item_id = job_id
                        params = {'min_cluster_size': min_size}
                        thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Perform Clustering", path_analysis, params))
                        st.session_state.active_thread = thread
                        thread.start()
                        st.info("Started HDBSCAN clustering.")
                    else: st.warning("Another pipeline is already running.")

        with st.expander("**3. Inspect a Data File**"):
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
                st_autorefresh(interval=2000, key="global_refresher")
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
        

        tab_plots, tab_logs, tab_status = st.tabs(["Visualisations", "Logs", "Resources"])
        
        with tab_plots:

            PLOTS = {
                "UMAP Projection": {"file": "snf_umap_projection.html", "header": "UMAP Projection of Relative Distance"},
                "HDBSCAN Clustering": {"file": "cluster_plot.html", "header": "HDBSCAN Clustering On Relative Distances"},
                "Cluster Stability": {"file": "stability_plot.html", "header": "Cluster Stability Analysis"},
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

        with tab_status:
            st_autorefresh(interval=1000, key="resource_refresher")
            colored_header(
                label="Resource Monitoring",
                description="Auto-refreshing snapshot of system resources.",
                color_name="green-70"
            )

            col1, col2, col3 = st.columns(3)
            
            cpu_usage = psutil.cpu_percent()
            col1.metric(label="System CPU Utilization", value=f"{cpu_usage}%")
            
            mem_info = psutil.virtual_memory()
            col2.metric(label="System RAM Usage", value=f"{(mem_info.total - mem_info.available) / (1024**3):.2f} GB")

            current_process = psutil.Process()
            process_mem_mb = current_process.memory_info().rss / (1024**2)
            col3.metric(label="App RAM Usage", value=f"{process_mem_mb:.2f} MB")
            
            style_metric_cards(border_left_color="#0068C9")