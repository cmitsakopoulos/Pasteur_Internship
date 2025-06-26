import streamlit as st
import pandas as pd
import time
import psutil
import os
import threading
from io import StringIO
import contextlib

from database_pipeline import (
    FileLogger, Pipeline, Concatenation, CDRComputation, 
    AntigenComputation, FlattenDuplicates, Write, RmPurificationTags, 
    AssignIDs, ComputeRelationships, WorkWithDatabase, CleanUp, 
    PreWalker, LevenshteinDistance, Walker
)
from function_dump import inspect_summary, inspect_verbose

RECIPES = {
    "Standard": [CleanUp, PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
    "Rerun": [PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
    "Distance Matrix": [PreWalker, Walker, Concatenation, LevenshteinDistance]
}

def build_pipeline(recipe: str, path: str, logger: FileLogger) -> Pipeline:
    logger.log(f"Building '{recipe}' pipeline for path: {path}")
    steps = []
    for step_class in RECIPES[recipe]:
        if step_class is PreWalker:
            steps.append(step_class(input_path=path, logger=logger))
        else:
            steps.append(step_class(logger=logger))
    return Pipeline(steps, logger)

def run_pipeline_worker(job_id, recipe, path):
    logger = None
    try:
        logger = FileLogger(job_id)
        pipeline = build_pipeline(recipe, path, logger)
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

def main():
    st.set_page_config(page_title="Pipeline Dashboard", layout="wide")

    if "text_output" not in st.session_state: st.session_state.text_output = ""
    if "running_job_id" not in st.session_state: st.session_state.running_job_id = None
    if "active_thread" not in st.session_state: st.session_state.active_thread = None

    control_col, dashboard_col = st.columns((1, 2))

    with control_col:
        st.title("🧬 Pipeline Control Panel")
        st.markdown("Configure and execute data processing and analysis tasks.")

        if st.session_state.running_job_id and st.session_state.active_thread:
            if not st.session_state.active_thread.is_alive():
                job_id = st.session_state.running_job_id
                log_filename = os.path.abspath(os.path.join("Internal_Files", f"{job_id}.log"))

                if os.path.exists(log_filename):
                    try:
                        with open(log_filename, 'r', encoding='utf-8') as f:
                            st.session_state.text_output = f.read()
                        os.remove(log_filename)
                        st.success("Pipeline run finished!")
                    except Exception as e:
                        st.error(f"Could not read or delete log file: {e}")
                else:
                    st.session_state.text_output = f"Job finished, but log file ('{log_filename}') was not found."

                st.session_state.running_job_id = None
                st.session_state.active_thread = None

        with st.expander("▶️ **1. Process Raw Data**", expanded=True):
            st.info("This pipeline processes raw experimental files into a structured format.")
            recipe = st.radio("Select Recipe", ('Standard', 'Rerun'),
                captions=["Clears prior results and runs fresh.", "Runs on new files only."],
                horizontal=True)
            path_run = st.text_input("Input Directory", value="./Internal_Files",
                help="Specify the directory containing the raw data to be processed.")

            if st.button("Execute Pipeline", type="primary", use_container_width=True):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id
                    abs_path = os.path.abspath(path_run)
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, recipe, abs_path))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info(f"Started '{recipe}' pipeline in the background.")
                else:
                    st.warning("Another pipeline is already running.")
        
        with st.expander("🧮 **2. Compute Distance Matrix**"):
            st.info("Calculates the Levenshtein distance between sequences in processed data.")
            available_dirs_dist = get_subdirectories()
            default_index_dist = available_dirs_dist.index("Internal_Files") if "Internal_Files" in available_dirs_dist else 0
            path_dist = st.selectbox("Select Source Directory", available_dirs_dist, index=default_index_dist, key="path_dist_select", help="Select the directory containing the processed CSV file.")
            if st.button("Calculate and Save Matrix", type="primary", use_container_width=True):
                if not st.session_state.running_job_id:
                    st.session_state.text_output = ""
                    job_id = f"job_{time.time()}"
                    st.session_state.running_job_id = job_id
                    abs_path_dist = os.path.abspath(path_dist)
                    thread = threading.Thread(target=run_pipeline_worker, args=(job_id, "Distance Matrix", abs_path_dist))
                    st.session_state.active_thread = thread
                    thread.start()
                    st.info("Started distance matrix calculation in the background.")
                else:
                    st.warning("Another pipeline is already running.")

        with st.expander("🔍 **3. Inspect a Data File**"):
            st.info("Quickly analyze the contents and structure of any CSV file.")
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
                        st.success("Quick summary generated.")
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
                        st.success("Detailed analysis complete.")
                    except Exception as e: st.error(f"Failed to run analysis: {e}")
            if not uploaded_file: st.caption("Upload a file to enable inspection options.")

    with dashboard_col:
        st.header("Activity & Results")
        tab_results, tab_status = st.tabs(["Results", "Resources"])
        with tab_results:
            st.subheader("Log / Text Output")
            if st.session_state.running_job_id:
                st.info("Pipeline is running in the background... The log will appear here upon completion.")
                with st.spinner("Processing..."):
                    while st.session_state.active_thread and st.session_state.active_thread.is_alive():
                        time.sleep(0.5)
                time.sleep(0.1)
                st.rerun()
            if st.session_state.text_output:
                st.code(st.session_state.text_output, language='bash')
            elif not st.session_state.running_job_id:
                st.info("Output from a pipeline run or file inspection will be shown here.")
        with tab_status:
            st.subheader("Live System Status")
            current_process = psutil.Process()
            if "monitoring" not in st.session_state: st.session_state.monitoring = False
            if st.button("Toggle Live Monitoring"):
                st.session_state.monitoring = not st.session_state.monitoring
                st.rerun() 
            if st.session_state.monitoring:
                st.success("Live monitoring is active.")
                placeholder = st.empty()
                while st.session_state.monitoring:
                    with placeholder.container():
                        st.markdown("##### System-Wide Resources")
                        col1, col2 = st.columns(2)
                        with col1:
                            cpu_usage = psutil.cpu_percent()
                            st.metric(label="System CPU Utilization", value=f"{cpu_usage}%")
                            st.progress(int(cpu_usage))
                        with col2:
                            mem_info = psutil.virtual_memory()
                            st.metric(label="System RAM Usage", value=f"{(mem_info.total - mem_info.available) / (1024**3):.2f} GB / {mem_info.total / (1024**3):.2f} GB")
                            st.progress(int(mem_info.percent))
                        st.markdown("##### Application-Specific Resources")
                        process_mem_mb = current_process.memory_info().rss / (1024**2)
                        st.metric(label="App RAM Usage", value=f"{process_mem_mb:.2f} MB")
                        st.caption("Represents the memory allocated to the Streamlit dashboard process (RSS).")
                    time.sleep(1) 
            else: 
                st.info("Live monitoring is inactive. Displaying a static resource snapshot.")
                st.markdown("##### System-Wide Resources")
                col1, col2 = st.columns(2)
                with col1:
                    cpu_usage = psutil.cpu_percent()
                    st.metric(label="System CPU Utilization", value=f"{cpu_usage}%")
                with col2:
                    mem_info = psutil.virtual_memory()
                    st.metric(label="System RAM Usage", value=f"{(mem_info.total - mem_info.available) / (1024**3):.2f} GB")
                st.markdown("##### Application-Specific Resources")
                process_mem_mb = current_process.memory_info().rss / (1024**2)
                st.metric(label="App RAM Usage", value=f"{process_mem_mb:.2f} MB")

if __name__ == "__main__":
    main()