# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 18:08:04 2025

@author: zengt
"""

import streamlit as st
import pandas as pd
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import os
import pickle
import plotly.express as px
import networkx as nx
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import io

script_dir = os.path.dirname(os.path.dirname(__file__))

st.set_page_config(page_title="TriNET", layout="centered",page_icon=":material/graph_5:")
st.markdown("### Functional protein complexes")
st.markdown("---")

# 读取数据
df_all = pickle.load(open(os.path.join(script_dir, "data/de_novo_TriNET_08_random.pkl"), 'rb'))
df_all = df_all.sort_values(by='Complex size',ascending=True)

# ====== 添加侧边栏统一阈值设置 ======
st.sidebar.markdown("## Filters")
score_threshold = st.sidebar.slider("Minimum interaction score", 0.8, 1.0, 0.8, 0.01)

df = df_all[df_all["TriNET score"] > score_threshold].copy()

st.markdown("#### TriNET results")

st.caption("📌 Tip: Click a row in the table below to view its TriNET details.")
st.info("Only a subset of interactions is shown below.")

# ---- 筛选列选择 ----
col_to_filter = st.selectbox("🔍 Select a column to filter", df.columns)

if pd.api.types.is_numeric_dtype(df[col_to_filter]):
    min_val = float(df[col_to_filter].min())
    max_val = float(df[col_to_filter].max())
    selected_range = st.slider(f"Filter {col_to_filter}", min_val, max_val, (min_val, max_val))
    filtered_df = df[(df[col_to_filter] >= selected_range[0]) & (df[col_to_filter] <= selected_range[1])]
else:
    input_val = st.text_input(f"Filter {col_to_filter}")
    if input_val:
        filtered_df = df[df[col_to_filter] == input_val]
    else:
        filtered_df = df.copy()

# ---- 构建 AgGrid ----
gb = GridOptionsBuilder.from_dataframe(filtered_df)
gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)  # 每页20行，可自定义
gb.configure_selection(selection_mode="single", use_checkbox=False)
gb.configure_grid_options(enableCellTextSelection=True)
grid_options = gb.build()

grid_response = AgGrid(
    filtered_df,
    gridOptions=grid_options,
    height=450,
    width='100%',
    update_mode=GridUpdateMode.MODEL_CHANGED,
    allow_unsafe_jscode=True,
    enable_enterprise_modules=False
)

# **新增：下载按钮**
csv = filtered_df.to_csv(index=False)
st.download_button(
    label="Download",
    icon=':material/download:',
    data=csv,
    file_name='TriNET.csv',
    mime='text/csv'
)

# ---- 获取选中行 ----
selected_rows = grid_response['selected_rows']

#%%

if isinstance(selected_rows, pd.DataFrame):
    selected_rows = selected_rows.to_dict(orient="records")

@st.cache_data
def load_tpi_comp_source():
    return pickle.load(open(os.path.join(script_dir, "data/tpi_comp_source.pkl"), 'rb'))

@st.cache_data
def load_de_novo_tpis():
    return pickle.load(open(os.path.join(script_dir, "data/de_novo_TPIs.pkl"), 'rb'))

if selected_rows:
    complex_genes_str = selected_rows[0].get("Complex", "")
    
    st.markdown("---")
    st.markdown("#### Network of {}".format(complex_genes_str))
    
    if complex_genes_str:
        genes = [g.strip() for g in complex_genes_str.split('-') if g.strip()]
        if not genes:
            st.warning("No valid genes found in selected interaction.")
        else:
            # 创建图
            G = nx.Graph()
    
            # 添加节点和边
            for g in genes:
                G.add_node(g, category='gene')
            for i in range(len(genes)):
                for j in range(i + 1, len(genes)):
                    G.add_edge(genes[i], genes[j], relation='gene-gene')
    
            # 设置节点/边样式
            node_color_map = {'gene': '#ffb6c1', 'drug': '#0a75ad'}
            node_shape_map = {'gene': 'o', 'drug': '*'}
            edge_color_map = {'gene-gene': '#ffb6c1', 'drug-gene': '#0a75ad'}
    
            pos = nx.kamada_kawai_layout(G)
    
            fig, ax = plt.subplots(figsize=(3, 3))
            edge_colors = [edge_color_map[G.edges[e].get('relation', 'gene-gene')] for e in G.edges()]
            nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=2, ax=ax)
    
            for node, (x, y) in pos.items():
                category = G.nodes[node].get('category', 'gene')
                shape = node_shape_map.get(category, 'o')
                color = node_color_map.get(category, '#cccccc')
                ax.scatter([x], [y], s=700, c=[color], marker=shape)
                ax.text(x, y + 0.3, node, fontsize=12, ha='center', va='center')
    
            ax.axis('off')
    
            buf = io.BytesIO()
            fig.savefig(buf, format="png", bbox_inches='tight')
            buf.seek(0)
            st.image(buf, width=300)
            plt.close(fig)
        
    # ---- 展示对应的新 dataframe ----
    st.markdown("---")
    st.markdown("#### Associated TPI profiles")
    
    ge_key = selected_rows[0].get("Data source", "")
    comp = complex_genes_str  # gene_list 就是 complex_genes_str.split("-")
    
    tpi_comp_source = load_tpi_comp_source()
    tpis = tpi_comp_source[ge_key][comp]
    
    # 读取数据
    df_all = load_de_novo_tpis()
    df_all = df_all.sort_values(by="TPCA-PCorr-TPI score", ascending=False)
    
    sub_df = df_all[df_all['TPI'].isin(tpis)].reset_index(drop=True)
    st.dataframe(sub_df)













