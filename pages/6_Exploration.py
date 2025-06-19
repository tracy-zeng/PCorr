# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 10:23:08 2025

@author: zengt
"""

import streamlit as st
import pandas as pd
import os
import pickle
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import re
import networkx as nx
import matplotlib.pyplot as plt
import io

script_dir = os.path.dirname(os.path.dirname(__file__))

st.set_page_config(page_title="Exploration", layout="centered", page_icon=":material/diagonal_line:")
st.markdown("### Hypothesis generating")
st.markdown("---")

query_type = st.radio("👇 Select query type:", ["Drug", "Gene"], horizontal=True, index=0)

def standardize_smi(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    clean_mol = rdMolStandardize.Cleanup(mol)
    return Chem.MolToSmiles(clean_mol, isomericSmiles=False, canonical=True)

@st.cache_data
def load_relation_data():
    df_ppi_all = pickle.load(open(os.path.join(script_dir, "data/de novo PPIs.pkl"), 'rb'))
    df_tpi_all = pickle.load(open(os.path.join(script_dir, "data/de novo TPIs.pkl"), 'rb'))
    df_trinet_all = pickle.load(open(os.path.join(script_dir, "data/de novo TriNET.pkl"), 'rb'))
    return df_ppi_all, df_tpi_all, df_trinet_all

@st.cache_data
def load_DGI_data():
    df_dgi_all = pickle.load(open(os.path.join(script_dir, "data/de novo DGIs.pkl"), 'rb'))
    return df_dgi_all

def draw_network_from_genes(gene_str, dgi_df):
    genes = [g.strip() for g in gene_str.split('-') if g.strip()]
    if not genes:
        st.warning("No valid genes found in selected interaction.")
        return

    G = nx.Graph()

    for g in genes:
        G.add_node(g, category='gene')

    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            G.add_edge(genes[i], genes[j], relation='gene-gene')

    dgi_sub = dgi_df[dgi_df['Gene'].isin(genes)]

    filtered_dgi_rows = []
    for gene in genes:
        dgi_gene = dgi_sub[dgi_sub['Gene'] == gene]
        if not dgi_gene.empty:
            dgi_top1 = dgi_gene.sort_values(by='PCorr-DGI score', ascending=False).head(1)
            if len(dgi_gene) > 1:
                st.info(f"{gene} interacts with {len(dgi_gene)} drug(s); showing top 1 by PCorr-DGI score.")
            filtered_dgi_rows.append(dgi_top1)
    if filtered_dgi_rows:
        dgi_sub_filtered = pd.concat(filtered_dgi_rows)
    else:
        dgi_sub_filtered = pd.DataFrame()

    for _, row in dgi_sub_filtered.iterrows():
        drug = row['Drug name'].split('; ')[0].capitalize()
        gene = row['Gene']
        G.add_node(drug, category='drug')
        G.add_edge(drug, gene, relation='drug-gene')

    pos = nx.kamada_kawai_layout(G)
    node_color_map = {'gene': '#ffb6c1', 'drug': '#0a75ad'}
    node_shape_map = {'gene': 'o', 'drug': '*'}
    edge_color_map = {'gene-gene': '#ffb6c1', 'drug-gene': '#0a75ad'}

    fig, ax = plt.subplots(figsize=(3, 3))
    edge_colors = [edge_color_map[G.edges[e]['relation']] for e in G.edges()]
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=2, ax=ax)
    for node, (x, y) in pos.items():
        category = G.nodes[node].get('category', 'gene')
        shape = node_shape_map.get(category, 'o')
        color = node_color_map.get(category, '#cccccc')
        ax.scatter([x], [y], s=700, c=[color], marker=shape)
        ax.text(x, y + 0.1, node, fontsize=12, ha='center', va='bottom')
    ax.axis('off')

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight')
    buf.seek(0)
    st.image(buf, width=300)
    plt.close(fig)

def display_relation_table(title, df, key, complex_col='Complex'):
    st.markdown(f"#### {title}")
    
    # 初始化 table_reload 为字典
    if 'table_reload' not in st.session_state:
        st.session_state['table_reload'] = {}
    if key not in st.session_state['table_reload']:
        st.session_state['table_reload'][key] = 0
    
    dynamic_key = f"{key}_{st.session_state['table_reload'][key]}"

    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_selection(selection_mode="single", use_checkbox=False)
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
    grid_options = gb.build()

    grid_response = AgGrid(
        df, 
        gridOptions=grid_options, 
        update_mode=GridUpdateMode.SELECTION_CHANGED, 
        height=300, 
        theme="streamlit", 
        key=dynamic_key
    )

    selected_rows = grid_response['selected_rows']
    if isinstance(selected_rows, pd.DataFrame):
        selected_rows = selected_rows.to_dict(orient="records")

    if selected_rows:
        gene_str = selected_rows[0][complex_col]
        draw_network_from_genes(gene_str, st.session_state.get("df_result", pd.DataFrame()))

        # 正确更新该表格 key 对应的刷新次数
        st.session_state['table_reload'][key] += 1

if query_type == "Drug":
    df_all = load_DGI_data()
    input_type = st.selectbox("Select input type:", ["Drug id", "Drug name", "SMILES"])
    user_query = st.text_input("Enter value:")

    if user_query:
        df_result = pd.DataFrame()
        if input_type == "Drug id":
            df_result = df_all[df_all["Drug id"] == user_query].reset_index(drop=True)
        elif input_type == "Drug name":
            if "name_normalized" not in df_all.columns:
                df_all["name_normalized"] = df_all["Drug name"].str.lower().apply(lambda x: re.sub(r'[^a-z0-9]', '', str(x)))
            user_query_lc = re.sub(r'[^a-z0-9]', '', user_query.lower())
            matched_df = df_all[df_all["name_normalized"].str.contains(user_query_lc, na=False)]
            if not matched_df.empty:
                drug_options = [f"{row['Drug id']} | {row['Drug name']}" for _, row in matched_df.iterrows()]
                selected_drugs = st.multiselect(f"{len(set(drug_options))} matched drugs:", sorted(set(drug_options)))
                selected_ids = [x.split(" | ")[0] for x in selected_drugs]
                df_selected = matched_df[matched_df["Drug id"].isin(selected_ids)]
                if not df_selected.empty:
                    df_result = df_selected.iloc[:, :-1].reset_index(drop=True)
        elif input_type == "SMILES":
            query_smi_std = standardize_smi(user_query)
            if query_smi_std:
                matched_df = df_all[df_all["SMILES"] == query_smi_std]
                if not matched_df.empty:
                    df_result = matched_df.reset_index(drop=True)

        if not df_result.empty:
            st.session_state["df_result"] = df_result
            df_ppi_all, df_tpi_all, df_trinet_all = load_relation_data()
            score_threshold = st.sidebar.slider("Minimum interaction score", 0.6, 1.0, 0.8, 0.01)

            def extract_genes_from_col(series):
                genes = set()
                for entry in series.dropna():
                    genes.update(entry.split('-'))
                return genes

            df_ppi_all = df_ppi_all[df_ppi_all["PCorr-PPI score"] > score_threshold]
            df_tpi_all = df_tpi_all[df_tpi_all["TPCA-PCorr-TPI score"] > score_threshold]
            df_trinet_all = df_trinet_all[df_trinet_all["TriNET score"] > score_threshold]
            df_result = df_result[df_result["PCorr-DGI score"] > score_threshold]

            genes_in_ppi = extract_genes_from_col(df_ppi_all["PPI"])
            genes_in_tpi = extract_genes_from_col(df_tpi_all["TPI"])
            genes_in_trinet = extract_genes_from_col(df_trinet_all["Complex"])
            genes_with_relations = genes_in_ppi.union(genes_in_tpi).union(genes_in_trinet)

            filter_related_only = st.checkbox("Filtering DGIs with PPI/TPI/TriNET interactions", value=False)
            df_result_display = df_result[df_result["Gene"].isin(genes_with_relations)].reset_index(drop=True) if filter_related_only else df_result.copy()

            gb = GridOptionsBuilder.from_dataframe(df_result_display)
            gb.configure_selection(selection_mode="single", use_checkbox=False)
            gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)
            grid_options = gb.build()

            st.markdown("#### Matched DGIs")
            grid_response = AgGrid(df_result_display, gridOptions=grid_options, update_mode=GridUpdateMode.SELECTION_CHANGED, height=300, theme="streamlit")
            selected_rows = grid_response['selected_rows']
            if isinstance(selected_rows, pd.DataFrame):
                selected_rows = selected_rows.to_dict(orient="records")
            if selected_rows:
                selected_gene = selected_rows[0].get("Gene")
                st.session_state['selected_drug'] = selected_rows[0].get("Drug name", "DemoDrug")
                if selected_gene:
                    st.markdown("---")
                    df_ppi_sub = df_ppi_all[df_ppi_all["PPI"].str.split('-').apply(lambda genes: selected_gene in genes)].rename(columns={"PPI": "Complex"})
                    
                    any_shown = False
                    if not df_ppi_sub.empty:
                        st.markdown(f"#### Interactions for {selected_gene}")
                        display_relation_table("PPIs", df_ppi_sub, key="ppi")
                        any_shown = True
                    df_tpi_sub = df_tpi_all[df_tpi_all["TPI"].str.split('-').apply(lambda genes: selected_gene in genes)].rename(columns={"TPI": "Complex"})
                    if not df_tpi_sub.empty:
                        st.markdown(f"#### Interactions for {selected_gene}")
                        display_relation_table("TPIs", df_tpi_sub, key="tpi")
                        any_shown = True
                    df_tri_sub = df_trinet_all[df_trinet_all["Complex"].str.split('-').apply(lambda genes: selected_gene in genes)]
                    if not df_tri_sub.empty:
                        st.markdown(f"#### Interactions for {selected_gene}")
                        display_relation_table("TriNET Complexes", df_tri_sub, key="trinet")
                        any_shown = True
                    if not any_shown:
                        st.info(f"No interactions for {selected_gene}.")
elif query_type == "Gene":
    # 加载数据
    df_ppi_all, df_tpi_all, df_trinet_all = load_relation_data()
    df_all = load_DGI_data()

    # 输入框
    gene_input = st.text_input("Enter a gene symbol:")

    if gene_input:
        # 设置阈值
        score_threshold = st.sidebar.slider("Minimum interaction score", 0.6, 1.0, 0.8, 0.01)
        gene = gene_input.strip()

        # 提取子集
        df_ppi_sub = df_ppi_all[
            (df_ppi_all["PCorr-PPI score"] > score_threshold) &
            (df_ppi_all["PPI"].str.split('-').apply(lambda genes: gene in genes))
        ].copy()

        df_tpi_sub = df_tpi_all[
            (df_tpi_all["TPCA-PCorr-TPI score"] > score_threshold) &
            (df_tpi_all["TPI"].str.split('-').apply(lambda genes: gene in genes))
        ].copy()

        df_tri_sub = df_trinet_all[
            (df_trinet_all["TriNET score"] > score_threshold) &
            (df_trinet_all["Complex"].str.split('-').apply(lambda genes: gene in genes))
        ].copy()

        # 可选 DGI 过滤
        filter_dgi_only = st.checkbox("Filtering interactions with DGIs", value=False)
        if filter_dgi_only:
            dgi_genes = set(df_all["Gene"])
            df_ppi_sub = df_ppi_sub[df_ppi_sub["PPI"].str.split('-').apply(lambda genes: any(g in dgi_genes for g in genes))]
            df_tpi_sub = df_tpi_sub[df_tpi_sub["TPI"].str.split('-').apply(lambda genes: any(g in dgi_genes for g in genes))]
            df_tri_sub = df_tri_sub[df_tri_sub["Complex"].str.split('-').apply(lambda genes: any(g in dgi_genes for g in genes))]

        # 存入 session_state，用于网络图展示
        st.session_state["df_result"] = df_all
        
        st.markdown("---")

        # 展示各表，先 rename 为 Complex
        any_shown = False
        if not df_ppi_sub.empty:
            df_ppi_sub = df_ppi_sub.rename(columns={"PPI": "Complex"})
            st.markdown(f"#### Interactions for {gene_input}")
            display_relation_table("PPIs", df_ppi_sub, key="gene_ppi")
            any_shown = True
        if not df_tpi_sub.empty:
            df_tpi_sub = df_tpi_sub.rename(columns={"TPI": "Complex"})
            st.markdown(f"#### Interactions for {gene_input}")
            display_relation_table("TPIs", df_tpi_sub, key="gene_tpi")
            any_shown = True
        if not df_tri_sub.empty:
            st.markdown(f"#### Interactions for {gene_input}")
            display_relation_table("TriNET Complexes", df_tri_sub, key="gene_trinet")
            any_shown = True
        if not any_shown:
            st.info(f"No interactions for {gene_input}.")






















