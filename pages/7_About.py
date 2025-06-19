# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 17:48:50 2025

@author: zengt
"""

import streamlit as st
import os

st.set_page_config(page_title="About",
                   layout="centered",
                   page_icon=":material/info:",)

st.title("About")

st.markdown("""Zeng et al. proposed PCorr, a modeling framework leveraging
            phenotypic correlation signals from large-scale cancer vulnerability
            profiles to predict functional protein–protein interactions (PCorr-PPI),
            ternary protein interactions (PCorr-TPI), and drug–target interactions (PCorr-DTI).
            Incorporating thermal proximity coaggregation (TPCA) signatures of protein assemblies,
            they developed an ensemble TPCA-PCorr-TPI classifier and further introduced the TriNET algorithm
            to enable de novo identification of functional protein complexes. They further developed
            L-PCorr by integrating lineage-specific correlation patterns to uncover context-dependent
            molecular interactions."""
            )

st.markdown("#### Organization:")
st.markdown("Southern University of Science and Technology")
st.markdown("#### Contributor:")
st.markdown("ZENG Ya (12250033@mail.sustech.edu.cn)")


st.markdown("\n")
script_dir = os.path.dirname(os.path.dirname(__file__))
st.image(os.path.join(script_dir, "images/logo.jpg"),width=300)






