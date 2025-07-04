import json
import os
import streamlit as st
import pandas as pd
from st_aggrid import AgGrid
from app_components.config import IMAGES_DIR, INTERNAL_FILES_DIR

def render_welcome_page():
    st.image(str(IMAGES_DIR / "logo.png"), width=230)
    st.title("Welcome")
    st.markdown("""
    This platform allows you to parse, process, and analyse Antibody-Antigen data, with interactive graphs.

    **To get started:**
    1.  Upload a sample `.csv` file of your own data below.
    2.  Map your raw file's columns to the required data fields; the pipeline uses a standardised set of columns to work with.
    3.  Save the mapping profile, this will then be used by the pipeline to understand your data.
    4.  Navigate to the **Analysis Dashboard**, you are now ready to do your analysis.
    """)

    st.header("Create read instructions")
    uploaded_file = st.file_uploader("Upload a sample CSV file of your dataset", type="csv")

    if uploaded_file:
        try:
            df_sample = pd.read_csv(uploaded_file)
            st.write("File Preview:")
            AgGrid(df_sample, height=250, fit_columns_on_grid_load=True)

            user_columns = list(df_sample.columns)
            st.subheader("Map your columns to the standardised fields:")

            column_map = {
                "pdb_name": st.selectbox("Column for PDB/Antibody ID:", user_columns),
                "antigen_sequence": st.selectbox("Column for Antigen Sequence:", user_columns),
                "heavy_cdr3": st.selectbox("Column for Heavy Chain CDR3:", user_columns),
                "light_cdr3": st.selectbox("Column for Light Chain CDR3:", user_columns),
                "dataset": st.selectbox("Column for Dataset/Origin (optional):", [None] + user_columns),
            }

            if st.button("Save Mapping Profile", type="primary", use_container_width=True):
                instructions_filepath = os.path.join(INTERNAL_FILES_DIR, "csv_instructions.json")
                
                with open(instructions_filepath, 'w') as f:
                    json.dump(column_map, f, indent=4)
                
                st.success(f"Mapping profile saved to `{instructions_filepath}`!")
                st.info("You can now proceed to the 'Analysis Dashboard' using the sidebar on the left.")

        except Exception as e:
            st.error(f"An error occurred: {e}")