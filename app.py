import streamlit as st
from streamlit_option_menu import option_menu

from app_components.welcome_page import render_welcome_page
from app_components.analysis_page import render_dashboard_page

def main():
    st.set_page_config(
        page_title="Antibody Analysis Platform",
        layout="wide"
    )

    with st.sidebar:
        selected = option_menu(
            "Main Menu", 
            ["Welcome", "Analysis Dashboard"], 
            icons=['house', 'bar-chart-line'], 
            menu_icon="cast", 
            default_index=0
        )

    if selected == "Welcome":
        render_welcome_page()
    elif selected == "Analysis Dashboard":
        render_dashboard_page()

if __name__ == "__main__":
    main()