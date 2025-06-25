import streamlit as st
import pandas as pd
import time
from io import StringIO
import contextlib
import psutil

from database_pipeline import (
    Pipeline, Walker, Concatenation, CDRComputation, AntigenComputation, FlattenDuplicates, Write, RmPurificationTags, AssignIDs, ComputeRelationships, WorkWithDatabase, CleanUp, PreWalker,)

from function_dump import (
    inspect_summary, inspect_verbose,)

RECIPES = {
    "Update Database": [PreWalker, Walker, WorkWithDatabase],
    "Standard": [PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
    "Rerun": [CleanUp, PreWalker, Walker, Concatenation, FlattenDuplicates, RmPurificationTags, CDRComputation, AntigenComputation,
        AssignIDs, ComputeRelationships, Write],
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
                st.session_state.latest_result_df = result_df

        with st.expander("🔍 **Inspect a File**"):
            uploaded_file = st.file_uploader("Upload a CSV for inspection", type="csv")
            if st.button("Inspect File", use_container_width=True):
                head_df, summary = inspect_dummy_file(uploaded_file, st)
                if head_df is not None:
                    st.session_state.latest_result_df = head_df
                    st.session_state.inspection_summary = summary
            
            st.button("Visualize Data (Dummy)", use_container_width=True, disabled=True)

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

            if "monitoring" not in st.session_state:
                st.session_state.monitoring = False

            if st.button("Toggle Live Monitoring"):
                st.session_state.monitoring = not st.session_state.monitoring

            if st.session_state.monitoring:
                st.success("Live monitoring is ON. The stats below will update every 2 seconds.")
            else:
                st.warning("Live monitoring is OFF. Click the button to start.")

            cpu_metric = st.empty()
            mem_metric = st.empty()
            cpu_chart = st.empty()
            mem_chart = st.empty()

            while st.session_state.monitoring:
                cpu_usage = psutil.cpu_percent(interval=1)
                mem_info = psutil.virtual_memory()
                mem_percent = mem_info.percent

                mem_used_gb = f"{mem_info.used / (1024**3):.2f}"
                mem_total_gb = f"{mem_info.total / (1024**3):.2f}"

                cpu_metric.metric(label="System CPU Utilization", value=f"{cpu_usage}%")
                mem_metric.metric(label="System RAM Usage (GB)", value=f"{mem_used_gb} / {mem_total_gb}")
                
                cpu_chart.progress(int(cpu_usage), text=f"CPU: {cpu_usage}%")
                mem_chart.progress(int(mem_percent), text=f"RAM: {mem_percent}%")
                time.sleep(1)
            else: 
                cpu_usage = psutil.cpu_percent()
                mem_info = psutil.virtual_memory()
                cpu_metric.metric(label="System CPU Utilization", value=f"{cpu_usage}%")
                mem_metric.metric(label="System RAM Usage (GB)", value=f"{mem_info.used / (1024**3):.2f} / {mem_info.total / (1024**3):.2f}")
                cpu_chart.progress(int(cpu_usage), text=f"CPU: {cpu_usage}%")
                mem_chart.progress(int(mem_info.percent), text=f"RAM: {mem_info.percent}%")


if __name__ == "__main__":
    main()