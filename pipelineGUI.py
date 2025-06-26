import streamlit as st
import pandas as pd
import time
from io import StringIO
import contextlib
import psutil

from database_pipeline import (
    Pipeline, Walker, Concatenation, CDRComputation, AntigenComputation, FlattenDuplicates, Write, RmPurificationTags, AssignIDs, ComputeRelationships, WorkWithDatabase, CleanUp, PreWalker, LevenshteinDistance)

from function_dump import (
    inspect_summary, inspect_verbose,)

RECIPES = {
    "Update Database": [PreWalker, Walker, WorkWithDatabase],
    "Standard": [PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
    "Rerun": [CleanUp, PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
    "Distance Matrix": [PreWalker, Walker, Concatenation, LevenshteinDistance]
}

def build_pipeline(recipe: str, path: str) -> Pipeline:
    print(f"Building '{recipe}' pipeline for path: {path}") 
    steps = [step(path) if step is PreWalker else step() for step in RECIPES[recipe]]
    return Pipeline(steps)

def run(recipe, path, log_area):
    """
    Runs the selected pipeline and captures all printed output.
    The final DataFrame is written to disk by the pipeline, not held in memory.
    """
    log_area.info(f"Starting '{recipe}' recipe on '{path}'...")
    log_capture_string = StringIO()
    with contextlib.redirect_stdout(log_capture_string):
        try:
            pipeline = build_pipeline(recipe, path)
            pipeline.run() # Run the pipeline, no DataFrame is returned
            print("\n--- PIPELINE COMPLETED ---")

        except Exception as e:
            print(f"\n--- A FATAL ERROR OCCURRED ---")
            print(str(e))
            log_area.error(f"An error occurred: {e}")

    full_log = log_capture_string.getvalue()
    log_area.success(f"Recipe '{recipe}' completed successfully!")
    return full_log

def main():
    st.set_page_config(
        page_title="Pipeline Dashboard",
        layout="wide"
    )

    # Simplified session state: only one variable for text output is needed.
    if "text_output" not in st.session_state:
        st.session_state.text_output = ""

    # Main Layout (2 columns) 
    control_col, dashboard_col = st.columns((1, 2))

    # Left Column: 
    with control_col:
        st.title("🧬 Pipeline Control Panel")
        st.markdown("Configure and execute data processing and analysis tasks.")

        with st.expander("▶️ **1. Process Raw Data**", expanded=True):
            st.info("This pipeline processes raw experimental files into a structured format.")
            recipe = st.radio(
                "Select Recipe",
                ('Standard', 'Rerun'),
                captions=[
                    "Run on new data.",
                    "Delete existing results and re-process all data."
                ],
                horizontal=True
            )
            path_run = st.text_input(
                "Input Directory",
                value="./Internal_Files",
                help="Specify the directory containing the raw data to be processed."
            )

            if st.button("Execute Pipeline", type="primary", use_container_width=True):
                st.session_state.text_output = ""
                log_text = run(recipe, path_run, st)
                st.session_state.text_output = log_text

        with st.expander("🧮 **2. Compute Distance Matrix**"):
            st.info("Calculates the Levenshtein distance between sequences in processed data.")
            path_dist = st.text_input(
                "Source Directory",
                value="./Internal_Files",
                key="dist_matrix_path",
                help="Provide the path to the directory containing the processed CSV file (e.g., 'concatenated.csv')."
            )
            if st.button("Calculate and Save Matrix", type="primary", use_container_width=True):
                st.session_state.text_output = ""
                log_text = run("Distance Matrix", path_dist, st)
                st.session_state.text_output = log_text
                st.success("Distance matrix calculation complete.")


        with st.expander("🔍 **3. Inspect a Data File**"):
            st.info("Quickly analyze the contents and structure of any CSV file.")
            uploaded_file = st.file_uploader(
                "Upload a data file for analysis",
                type=["csv"],
                help="Upload a processed or external CSV file to generate a high-level summary or a detailed column-by-column analysis."
            )

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
                    except Exception as e:
                        st.error(f"Failed to generate summary: {e}")

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
                    except Exception as e:
                        st.error(f"Failed to run analysis: {e}")

            if not uploaded_file:
                st.caption("Upload a file to enable inspection options.")

    # Right Column:
    with dashboard_col:
        st.header("Activity & Results")

        tab_results, tab_status = st.tabs(["Results", "Resources"])

        with tab_results:
            st.subheader("Log / Text Output")
            
            # This is now the ONLY output area.
            if st.session_state.text_output:
                 st.code(st.session_state.text_output, language='bash')
            else:
                st.info("Output from a pipeline run or file inspection will be shown here.")

        with tab_status:
            st.subheader("Live System Status")
            # ... The rest of the tab_status code is unchanged ...
            current_process = psutil.Process()
            if "monitoring" not in st.session_state:
                st.session_state.monitoring = False
            if st.button("Toggle Live Monitoring"):
                st.session_state.monitoring = not st.session_state.monitoring
                st.rerun() 
            if st.session_state.monitoring:
                st.success("Live monitoring is active. System metrics will refresh periodically.")
            else:
                st.info("Live monitoring is inactive. Displaying a static resource snapshot. Activate monitoring for real-time updates.")
            st.markdown("##### System-Wide Resources")
            col1, col2 = st.columns(2)
            with col1:
                cpu_metric = st.empty()
                cpu_chart = st.empty()
            with col2:
                mem_metric = st.empty()
                mem_chart = st.empty()
            st.markdown("##### Application-Specific Resources")
            app_mem_metric = st.empty()
            st.caption("Represents the memory allocated to the Streamlit dashboard process (Resident Set Size).")
            if st.session_state.monitoring:
                while True:
                    cpu_usage = psutil.cpu_percent(interval=1)
                    mem_info = psutil.virtual_memory()
                    process_mem_mb = current_process.memory_info().rss / (1024**2)
                    cpu_metric.metric(label="System CPU Utilization", value=f"{cpu_usage}%")
                    cpu_chart.progress(int(cpu_usage))
                    mem_metric.metric(label="System RAM Usage", value=f"{(mem_info.total - mem_info.available) / (1024**3):.2f} GB / {mem_info.total / (1024**3):.2f} GB")
                    mem_chart.progress(int(mem_info.percent))
                    app_mem_metric.metric(label="App RAM Usage", value=f"{process_mem_mb:.2f} MB")
                    time.sleep(1) 
            else: 
                cpu_usage = psutil.cpu_percent()
                mem_info = psutil.virtual_memory()
                process_mem_mb = current_process.memory_info().rss / (1024**2)
                cpu_metric.metric(label="System CPU Utilization", value=f"{cpu_usage}%")
                cpu_chart.progress(int(cpu_usage))
                mem_metric.metric(label="System RAM Usage", value=f"{(mem_info.total - mem_info.available) / (1024**3):.2f} GB / {mem_info.total / (1024**3):.2f} GB")
                mem_chart.progress(int(mem_info.percent))
                app_mem_metric.metric(label="App RAM Usage", value=f"{process_mem_mb:.2f} MB")

if __name__ == "__main__":
    main()