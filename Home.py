# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 14:16:09 2025

@author: zengt
"""

import streamlit as st
from streamlit.logger import get_logger
import os

LOGGER = get_logger(__name__)


def run():
    st.set_page_config(
        page_title="PCorr",
        page_icon=":material/forest:",
        layout="centered",
        initial_sidebar_state="auto",
    )
    
    script_dir = os.path.dirname(__file__)

    # 创建两列，第一列放内容，第二列放 logo 图片
    col1, col2 = st.columns([3, 1])  # 比例可调整为更靠右：如 [4,1]

    with col2:
        st.image(os.path.join(script_dir, "images/logo.jpg"), width=300)  # 缩小一点更美观

    with col1:
        st.markdown("<h1 style='font-weight:bold;'>PCorr</h1>", unsafe_allow_html=True)
        st.markdown("---")
        st.markdown(
            """
            ### Discovering molecular interactions from diverse perturbation data.
            """
        )

    st.markdown(" ")

    st.image(os.path.join(script_dir, "images/Graphical abstract.svg"), width=500)

if __name__ == "__main__":
    run()








