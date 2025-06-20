# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 13:51:09 2025

@author: zengt
"""

import streamlit as st
import pandas as pd
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import os
import pickle
import plotly.express as px

script_dir = os.path.dirname(os.path.dirname(__file__))

st.set_page_config(page_title="PCorr-PPI", layout="centered",page_icon=":material/diagonal_line:")
st.markdown("### Functional protein-protein interactions (PPIs)")
st.markdown("---")

# 读取数据
df_all = pickle.load(open(os.path.join(script_dir, "data/de_novo_PPIs_08.pkl"), 'rb'))
df_all = df_all.sort_values(by="PCorr-PPI score", ascending=False)

# ====== 添加侧边栏统一阈值设置 ======
st.sidebar.markdown("## Filters")
score_threshold = st.sidebar.slider("Minimum interaction score", 0.6, 1.0, 0.8, 0.01)

df = df_all[df_all["PCorr-PPI score"] > score_threshold].copy()

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

st.markdown("#### PCorr-PPI results")

st.caption("📌 Tip: Click a row in the table below to view its PCorr signal plot.")

# ---- 构建 AgGrid ----
gb = GridOptionsBuilder.from_dataframe(filtered_df)
gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)  # 每页20行，可自定义
gb.configure_selection(selection_mode="single", use_checkbox=False)
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
    file_name='pcorr_PPIs.csv',
    mime='text/csv'
)

# ---- 获取选中行 ----
selected_rows = grid_response['selected_rows']

#%%

ge_matrix = pickle.load(open(os.path.join(script_dir, "data/ge_matrix.pkl"), 'rb'))
ccl_lineage = pickle.load(open(os.path.join(script_dir, "data/ccl_lineage.pkl"), 'rb'))
lineage_ccl = pickle.load(open(os.path.join(script_dir, "data/lineage_ccl.pkl"), 'rb'))
lineage_color = pickle.load(open(os.path.join(script_dir, "data/lineage_color.pkl"), 'rb'))

ge_cut = {}
ge_cut['CRISPR'] = -0.237
ge_cut['RNAi'] = -0.325

def list_count(lst):
    count_dict = {}
    set_ = list(set(lst))
    for i in set_:
        count = lst.count(i)
        count_dict.update({i: count})
    return count_dict

def correlation_plot_plotly(g1, g2, ge_key, highlight_by=None):
    ge_c = ge_cut[ge_key]
    ge_1 = ge_matrix[ge_key][g1].dropna()
    ge_2 = ge_matrix[ge_key][g2].dropna()
    ccls = list(set(ge_1.index) & set(ge_2.index))

    # 构建数据
    data = []
    for ccl in ccls:
        lineage = ccl_lineage.get(ccl, "Unknown")
        data.append({
            "CCL": ccl,
            "Lineage": lineage,
            "X": ge_1[ccl],
            "Y": ge_2[ccl]
        })

    df_plot = pd.DataFrame(data)

    use_highlight = False
    if highlight_by:
        if highlight_by in df_plot["CCL"].values:
            df_plot["State"] = df_plot["CCL"].apply(lambda x: "Highlight" if x == highlight_by else "Other")
            use_highlight = True
        elif highlight_by in df_plot["Lineage"].values:
            df_plot["State"] = df_plot["Lineage"].apply(lambda x: "Highlight" if x == highlight_by else "Other")
            use_highlight = True
        else:
            st.warning(f'"{highlight_by}" not found in CCL or Lineage.')
            df_plot["State"] = "Highlight"
            use_highlight = False
    else:
        df_plot["State"] = "Highlight"  # 全部按 lineage 上色
        use_highlight = False

    # 分配颜色
    if use_highlight:
        df_plot["Color"] = df_plot.apply(
            lambda row: lineage_color[row["Lineage"]] if row["State"] == "Highlight" else "lightgray", axis=1)
        color_column = "State"
        color_values = df_plot["State"]
        color_map = {
            "Highlight": "red",  # 实际不生效，但必须占位
            "Other": "lightgray"
        }
        legend_title = "State"
    else:
        df_plot["Color"] = df_plot["Lineage"].map(lineage_color)
        color_column = "Lineage"
        color_values = df_plot["Lineage"]
        color_map = lineage_color
        legend_title = "Lineage"

    # 绘图
    fig = px.scatter(
        df_plot,
        x="X",
        y="Y",
        color=color_values,
        color_discrete_map=color_map,
        hover_data={"CCL": True, "Lineage": True, "State": use_highlight},
        labels={
            "X": f"Gene effect of {g1} ({ge_key})",
            "Y": f"Gene effect of {g2} ({ge_key})",
            color_column: legend_title
        }
    )

    fig.update_traces(marker=dict(size=13, line=dict(width=0, color='gray'), opacity=0.8))
    fig.add_vline(x=ge_c, line_dash="dash", line_color="black")
    fig.add_hline(y=ge_c, line_dash="dash", line_color="black")

    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    fig.update_layout(
        title='',
        width=800,
        height=600,
        legend_title_text=legend_title,
        legend=dict(
            x=1.05,
            y=1,
            bgcolor='rgba(0,0,0,0)',
            borderwidth=0
        ),
        margin=dict(l=40, r=160, t=60, b=40),
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True),
        title_font_size=20
    )

    return fig

if isinstance(selected_rows, pd.DataFrame):
    selected_rows = selected_rows.to_dict(orient="records")

has_selection = isinstance(selected_rows, list) and len(selected_rows) > 0

if not has_selection:
    # 显示示例图
    example_ge_key = 'CRISPR'
    gene_a, gene_b = 'NFKB1','RELA'
    
    st.markdown("---")
    st.markdown("#### Example PCorr signal of {}-{} ({})".format(gene_a,gene_b,example_ge_key))

    highlight_term = st.text_input("🔍 Search for a lineage or cancer cell line (CCL) to highlight", value="", placeholder="e.g. Lymphoid or KMM1")

    fig = correlation_plot_plotly(gene_a, gene_b, example_ge_key, highlight_by=highlight_term.strip() if highlight_term else None)
    st.plotly_chart(fig, use_container_width=False)
    
else:
    # 选中时显示对应图
    selected_row = selected_rows[0]
    ppi_str = selected_row.get("PPI", "")
    ge_key = selected_row.get("Data source", "")

    gene_a, gene_b = ppi_str.split("-")
    
    if ge_key == 'CRISPR; RNAi':
        st.markdown("---")
        st.markdown("#### PCorr signal of {}-{} (CRISPR)".format(gene_a,gene_b))
        
        highlight_c = st.text_input("🔍 Search for a lineage or cancer cell line (CCL) to highlight", value="", placeholder="e.g. Lymphoid or KMM1",
                                    key='highlight_c')
        
        fig = correlation_plot_plotly(gene_a, gene_b, 'CRISPR', highlight_by=highlight_c.strip() if highlight_c else None)
        st.plotly_chart(fig, use_container_width=False)
        
        
        st.markdown("---")
        st.markdown("#### PCorr signal of {}-{} (RNAi)".format(gene_a,gene_b))
        
        highlight_r = st.text_input("🔍 Search for a lineage or cancer cell line (CCL) to highlight", value="", placeholder="e.g. Lymphoid or KMM1",
                                    key='highlight_r')
        
        fig = correlation_plot_plotly(gene_a, gene_b, 'RNAi', highlight_by=highlight_r.strip() if highlight_r else None)
        st.plotly_chart(fig, use_container_width=False)
    else:
        st.markdown("---")
        st.markdown("#### PCorr signal of {}-{} ({})".format(gene_a,gene_b,ge_key))
        
        highlight_term = st.text_input("🔍 Search for a lineage or cancer cell line (CCL) to highlight", value="", placeholder="e.g. Lymphoid or KMM1")
        
        fig = correlation_plot_plotly(gene_a, gene_b, ge_key, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
























