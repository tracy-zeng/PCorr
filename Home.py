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
        layout="wide",
        initial_sidebar_state="expanded",
    )
    script_dir = os.path.dirname(__file__)
    st.write("# PCorr exploration starts!")

    st.markdown(
        """
        PCorr——a modeling framework leveraging phenotypic correlation signals 
        from large-scale cancer vulnerability profiles to predict functional interactions.
        """
        )
    st.markdown(" ")
    image_path = os.path.join(script_dir, "images/pcorr.svg")
    st.image(image_path,width=300)

if __name__ == "__main__":
    run()








