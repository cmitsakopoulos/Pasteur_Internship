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
    log_area.info(f"Starting '{recipe}' recipe on '{path}'...")
    log_capture_string = StringIO()
    with contextlib.redirect_stdout(log_capture_string):
        try:
            pipeline = build_pipeline(recipe, path)
            just_run_the_pipe = pipeline.run() 
            print("\n--- PIPELINE COMPLETED ---")

        except Exception as e:
            print(f"\n--- A FATAL ERROR OCCURRED ---")
            print(str(e))
            log_area.error(f"An error occurred: {e}")
            return None

    full_log = log_capture_string.getvalue()
    st.session_state.log_output = full_log.splitlines() 
    log_area.success(f"Recipe '{recipe}' completed successfully!")
    return just_run_the_pipe

def inspect_dummy_file(uploaded_file, log_area):
    """Placeholder for the 'inspect' command."""
    if not uploaded_file:
        log_area.warning("⚠️ Please upload a file to inspect.")
        return None, None
        
    log_area.info(f"🔍 Inspecting file: {uploaded_file.name}...")
    with st.spinner("Analyzing CSV structure..."):
        time.sleep(1.5)
        try:
            df = pd.read_csv(uploaded_file)
            summary = f"""
File: {uploaded_file.name}
Size: {uploaded_file.size} bytes
Rows: {len(df)}
Columns: {len(df.columns)}
---
Column Names:
{', '.join(df.columns)}
            """
            log_area.success("Inspection complete.")
            return df.head(), summary
        except Exception as e:
            log_area.error(f"Failed to read or inspect file: {e}")
            return None, None

def main():
    st.set_page_config(
        page_title="Pipeline Dashboard",
        layout="wide"
    )

    if "log_output" not in st.session_state:
        st.session_state.log_output = []
    if "latest_result_df" not in st.session_state:
        st.session_state.latest_result_df = pd.DataFrame()
    if "inspection_summary" not in st.session_state:
        st.session_state.inspection_summary = ""


    # Main Layout (2 columns) 
    control_col, dashboard_col = st.columns((1, 2))

    # Left Column: 
    with control_col:
        st.title("🧬 Pipeline Control Panel")
        st.markdown("Configure and execute pipeline tasks.")
        with st.expander("▶️ **Run a Recipe**", expanded=True):
            recipe = st.radio(
                "Select a preconfigured Recipe",
                ('Standard', 'Rerun'),
                captions=["Standard run.", "Wipe results and re-run."],
                horizontal=True
            )
            path_run = st.text_input("Input Directory", value="./Internal_Files", help="Directory containing raw data.")
            
            if st.button("Execute Pipeline", type="primary", use_container_width=True):
                st.session_state.log_output = []
                def log_to_state(message):
                    st.session_state.log_output.append(message)

                result_df = run(recipe, path_run, st)
                if result_df is not None:
                    st.session_state.latest_result_df = result_df

        with st.expander("🧮 **Compute Distance Matrix**"):
            path_dist = st.text_input("Source Directory", value="./Internal_Files", key="dist_matrix_path", help="")
            if st.button("Calculate and Save Matrix", type="primary", use_container_width=True):
                st.session_state.log_output = []
                st.session_state.latest_result_df = pd.DataFrame()
                result = run("Distance Matrix", path_dist, st)
                if result is not None:
                    st.session_state.latest_result_df = result
                    st.success("Distance matrix calculation complete. Check logs and results.")
                else:
                    st.error("Distance matrix calculation failed. Please check the logs for errors.")

        with st.expander("🔍 **Inspect a File**", expanded=True):
            uploaded_file = st.file_uploader(
                "Upload a data file for analysis",
                type=["csv"],
                help="Upload a CSV file to generate either a high-level summary or a detailed column-by-column analysis."
            )

            col1, col2 = st.columns(2)

            with col1:
                if st.button("Quick Summary", use_container_width=True, disabled=not uploaded_file):
                    st.session_state.latest_result_df = pd.DataFrame()
                    st.session_state.inspection_summary = ""
                    
                    log_capture_string = StringIO()
                    try:
                        with contextlib.redirect_stdout(log_capture_string):
                            df = pd.read_csv(uploaded_file)
                            inspect_summary(df) 
                        
                        st.session_state.inspection_summary = log_capture_string.getvalue()
                        st.success("Quick summary generated.")
                    except Exception as e:
                        st.error(f"Failed to generate summary: {e}")

            with col2:
                if st.button("Detailed Analysis", use_container_width=True, disabled=not uploaded_file):
                    st.session_state.latest_result_df = pd.DataFrame()
                    st.session_state.inspection_summary = ""
                    
                    log_capture_string = StringIO()
                    try:
                        with contextlib.redirect_stdout(log_capture_string):
                            df = pd.read_csv(uploaded_file)
                            inspect_verbose(df) 
                        
                        st.session_state.inspection_summary = log_capture_string.getvalue()
                        st.success("Detailed analysis complete.")
                    except Exception as e:
                        st.error(f"Failed to run analysis: {e}")

            if not uploaded_file:
                st.info("Please upload a file to enable inspection.")

    # Right Column:
    with dashboard_col:
        st.header("Activity & Results")

        tab_log, tab_results, tab_status = st.tabs(["Output Log", "Results", "Resources"])

        with tab_log:
            st.subheader("Real-time Log")
            log_container = st.container(height=400)
            if st.session_state.log_output:
                for line in st.session_state.log_output:
                    log_container.text(line)
            else:
                 log_container.info("Updates on pipeline operations will appear here...")

        with tab_results:
            st.subheader("Latest Data Output")
            if not st.session_state.latest_result_df.empty:
                st.dataframe(st.session_state.latest_result_df, use_container_width=True)
            else:
                st.info("Results from a pipeline run or file inspection will be shown here.")
            
            if st.session_state.inspection_summary:
                 st.subheader("Inspection Summary")
                 st.code(st.session_state.inspection_summary, language='bash')

        with tab_status:
            st.subheader("Live System Status")

            # Get the current python process to track its own memory usage
            current_process = psutil.Process()

            if "monitoring" not in st.session_state:
                st.session_state.monitoring = False

            if st.button("Toggle Live Monitoring"):
                st.session_state.monitoring = not st.session_state.monitoring
                st.rerun() 

            # More professional user messages
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
            # More professional caption
            st.caption("Represents the memory allocated to the Streamlit dashboard process (Resident Set Size).")

            if st.session_state.monitoring:
                while True:
                    # System-wide stats
                    cpu_usage = psutil.cpu_percent(interval=1)
                    mem_info = psutil.virtual_memory()

                    # Process-specific RAM (RSS: Resident Set Size)
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