# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 18:04:06 2025

@author: zengt
"""

import streamlit as st
import pandas as pd
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import os
import pickle
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw
import io

script_dir = os.path.dirname(os.path.dirname(__file__))

st.set_page_config(page_title="L-PCorr", layout="centered",page_icon=":material/my_location:")
st.markdown("### Lineage-specific PCorr signals")
st.markdown("---")

# 设置 interaction_type 的选项及默认值
interaction_options = ["PPI", "TPI", "DGI"]
default_index = interaction_options.index("PPI")

# interaction type 单选按钮
interaction_type = st.radio("👇 Select an interaction type:", options=interaction_options, index=default_index, horizontal=True)

# 根据选择读取不同的数据
if interaction_type == "TPI":
    df_all = pickle.load(open(os.path.join(script_dir, "data/ls_tpi_08.pkl"), 'rb'))
    df_all = df_all.drop_duplicates().reset_index(drop=True)
    df_all = df_all.sort_values(by="TPCA-PCorr-TPI score", ascending=False)
else:
    df_all = pickle.load(open(os.path.join(script_dir, "data/ls_{}_08.pkl".format(interaction_type.lower())), 'rb'))
    df_all = df_all.drop_duplicates().reset_index(drop=True)
    df_all = df_all.sort_values(by="PCorr-{} score".format(interaction_type), ascending=False)

# ====== 添加侧边栏统一阈值设置 ======
st.sidebar.markdown("## Filters")
score_threshold = st.sidebar.slider("Minimum interaction score", 0.8, 1.0, 0.8, 0.01)

if interaction_type == "TPI":
    df = df_all[df_all["TPCA-PCorr-TPI score"] > score_threshold].copy()
else:
    df = df_all[df_all["PCorr-{} score".format(interaction_type)] > score_threshold].copy()

#%%
# ---- 筛选列选择 ----
col_to_filter = st.selectbox("🔍 Select a column to filter", df.columns)

if pd.api.types.is_numeric_dtype(df[col_to_filter]):
    min_val = float(df[col_to_filter].min())
    max_val = float(df[col_to_filter].max())
    selected_range = st.slider(f"Filter {col_to_filter}", min_val, max_val, (min_val, max_val))
    filtered_df = df[(df[col_to_filter] >= selected_range[0]) & (df[col_to_filter] <= selected_range[1])]
else:
    input_val = st.text_input(f"Filter {col_to_filter}", key=f"text_input_all_{col_to_filter}")
    if input_val:
        filtered_df = df[df[col_to_filter] == input_val]
    else:
        filtered_df = df.copy()

st.markdown("#### All L-PCorr results")

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
    file_name='L-PCorr_{}s.csv'.format(interaction_type),
    mime='text/csv'
)

#%%

df_l = pickle.load(open(os.path.join(script_dir, "data/ls_{}_bylineage_08.pkl".format(interaction_type.lower())), 'rb'))

# lineage 选择框（默认选择 Skin）
lineages = sorted(list(df_l.keys()))
default_index = lineages.index('Myeloid')

l = st.selectbox("Select a lineage:", lineages, index=default_index)

df_all = df_l[l].iloc[:,1:]
df_all = df_all.drop_duplicates().reset_index(drop=True)

# 默认筛选 > 0.8
if interaction_type == "TPI":
    df = df_all[df_all["TPCA-PCorr-TPI score"] > score_threshold].copy()
else:
    df = df_all[df_all["PCorr-{} score".format(interaction_type)] > score_threshold].copy()

#%%
# ---- 筛选列选择 ----
col_to_filter = st.selectbox("🔍 Select a column to filter", df.columns)

if pd.api.types.is_numeric_dtype(df[col_to_filter]):
    min_val = float(df[col_to_filter].min())
    max_val = float(df[col_to_filter].max())
    selected_range = st.slider(f"Filter {col_to_filter}", min_val, max_val, (min_val, max_val))
    filtered_df = df[(df[col_to_filter] >= selected_range[0]) & (df[col_to_filter] <= selected_range[1])]
else:
    input_val = st.text_input(f"Filter {col_to_filter}", key=f"text_input_lineage_{col_to_filter}")
    if input_val:
        filtered_df = df[df[col_to_filter] == input_val]
    else:
        filtered_df = df.copy()

st.markdown("### {}-PCorr results".format(l))
st.caption("📌 Tip: Click a row in the table below to view its {}-PCorr signal plot(s).".format(l))

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
    file_name='{}-PCorr_{}s.csv'.format(l,interaction_type),
    mime='text/csv'
)

# ---- 获取选中行 ----
selected_rows = grid_response['selected_rows']

#%%

dr_matrix = pickle.load(open(os.path.join(script_dir, "data/dr_matrix.pkl"), 'rb'))
ge_matrix = pickle.load(open(os.path.join(script_dir, "data/ge_matrix.pkl"), 'rb'))
ccl_lineage = pickle.load(open(os.path.join(script_dir, "data/ccl_lineage.pkl"), 'rb'))
lineage_ccl = pickle.load(open(os.path.join(script_dir, "data/lineage_ccl.pkl"), 'rb'))
lineage_color = pickle.load(open(os.path.join(script_dir, "data/lineage_color.pkl"), 'rb'))

id_name = pickle.load(open(os.path.join(script_dir, "data/id_name.pkl"), 'rb'))
name_id = pickle.load(open(os.path.join(script_dir, "data/name_id.pkl"), 'rb'))

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

def PPI_plotly(g1, g2, ge_key, lineage_filter, highlight_by=None):
    ge_c = ge_cut[ge_key]
    ge_1 = ge_matrix[ge_key][g1].dropna()
    ge_2 = ge_matrix[ge_key][g2].dropna()
    ccls = list(set(ge_1.index) & set(ge_2.index))

    # 仅保留指定 lineage 的 CCLs
    ccls = [ccl for ccl in ccls if ccl_lineage.get(ccl) == lineage_filter]

    data = []
    for ccl in ccls:
        data.append({
            "CCL": ccl,
            "X": ge_1[ccl],
            "Y": ge_2[ccl],
            "Lineage": ccl_lineage.get(ccl)
        })
    df_plot = pd.DataFrame(data)

    # 检查 highlight_by 是否在有效 CCL 中
    use_highlight = False
    if highlight_by:
        if highlight_by in df_plot["CCL"].values:
            df_plot["State"] = df_plot["CCL"].apply(lambda x: "Highlight" if x == highlight_by else "Other")
            use_highlight = True
        else:
            import streamlit as st
            st.warning(f'"{highlight_by}" not found in CCLs of {lineage_filter}.')
            highlight_by = None
            use_highlight = False

    if use_highlight:
        # Highlight 特定 CCL，其他为灰色
        df_plot["Color"] = df_plot["State"].apply(lambda x: "red" if x == "Highlight" else "lightgray")
        color_column = "State"
        color_map = {"Highlight": "red", "Other": "lightgray"}
        legend_title = "State"
    else:
        # 没有高亮，按 lineage_color 上色
        df_plot["Color"] = df_plot["Lineage"].map(lambda x: lineage_color.get(x, "gray"))
        color_column = "Lineage"
        color_map = lineage_color
        legend_title = "Lineage"

    fig = px.scatter(
        df_plot,
        x="X",
        y="Y",
        color=color_column,
        color_discrete_map=color_map,
        hover_data=["CCL", "Lineage"],
        labels={
            "X": f"Gene effect of {g1} ({ge_key})",
            "Y": f"Gene effect of {g2} ({ge_key})"
        },
        category_orders={color_column: sorted(df_plot[color_column].unique())}
    )

    fig.update_traces(marker=dict(size=13, line=dict(width=0, color='gray'), opacity=0.8))
    fig.add_vline(x=ge_c, line_dash="dash", line_color="black")
    fig.add_hline(y=ge_c, line_dash="dash", line_color="black")
    fig.update_yaxes(scaleanchor="x", scaleratio=1)

    fig.update_layout(
        width=800,
        height=600,
        margin=dict(l=40, r=160, t=60, b=40),
        title='{}-{}'.format(g1,g2),
        legend_title_text=legend_title
    )

    return fig

def DGI_plotly(drug, gene, dr_key, ge_key, lineage_filter, highlight_by=None):
    dr = dr_matrix[dr_key][drug].dropna()
    
    ge_c = ge_cut[ge_key]
    ge = ge_matrix[ge_key][gene].dropna()
    
    ccls = list(set(dr.index) & set(ge.index))
    
    # 仅保留指定 lineage 的 CCLs
    ccls = [ccl for ccl in ccls if ccl_lineage.get(ccl) == lineage_filter]
    
    # 构建数据
    data = []
    for ccl in ccls:
        lineage = ccl_lineage.get(ccl, "Unknown")
        data.append({
            "CCL": ccl,
            "Lineage": lineage,
            "X": dr[ccl],
            "Y": ge[ccl]
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
            "X": "Drug response of {} ({})".format(sorted(list(id_name[drug]))[0],dr_key),
            "Y": f"Gene effect of {gene} ({ge_key})",
            color_column: legend_title
        }
    )

    fig.update_traces(marker=dict(size=13, line=dict(width=0, color='gray'), opacity=0.8))
    fig.add_hline(y=ge_c, line_dash="dash", line_color="black")
    
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
        xaxis=dict(showgrid=True),  # 自动适应 x 轴范围
        yaxis=dict(showgrid=True),
        title_font_size=20
    )
    return fig

def draw_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            img = Draw.MolToImage(mol, size=(300, 300), kekulize=True)
            buf = io.BytesIO()
            img.save(buf, format='PNG')
            buf.seek(0)
            st.image(buf, caption='', width=300)
        else:
            st.warning(f"Invalid SMILES: {smiles}")
    except Exception as e:
        st.error(f"Error generating structure: {e}")

if isinstance(selected_rows, pd.DataFrame):
    selected_rows = selected_rows.to_dict(orient="records")

has_selection = isinstance(selected_rows, list) and len(selected_rows) > 0

if interaction_type == 'TPI':
    if not has_selection:
        # 显示示例图
        gene_1,gene_2,gene_3,example_ge_key="ABL1","BCR",'GRB2','CRISPR'
        example_l = 'Myeloid'
        
        st.markdown("---")
        st.markdown("#### Example {}-PCorr signal of {}-{}-{} ({})".format(example_l,gene_1,gene_2,gene_3,example_ge_key))
    
        highlight_term = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="")
        
        fig = PPI_plotly(gene_1, gene_2, example_ge_key, lineage_filter=example_l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
        
        fig = PPI_plotly(gene_1, gene_3, example_ge_key, lineage_filter=example_l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
        
        fig = PPI_plotly(gene_2, gene_3, example_ge_key, lineage_filter=example_l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
    else:
        # 选中时显示对应图
        selected_row = selected_rows[0]
        tpi_str = selected_row.get(interaction_type, "")
        ge_key = selected_row.get("Data source", "")
    
        gene_a, gene_b, gene_c = tpi_str.split("-")
        
        st.markdown("---")
        st.markdown("#### {}-PCorr signal of {}-{}-{} ({})".format(l,gene_a,gene_b,gene_c,ge_key))
    
        highlight_term = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="")
        
        fig = PPI_plotly(gene_a, gene_b, ge_key, lineage_filter=l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
        
        fig = PPI_plotly(gene_a, gene_c, ge_key, lineage_filter=l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
        
        fig = PPI_plotly(gene_b, gene_c, ge_key, lineage_filter=l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
elif interaction_type == 'PPI':
    if not has_selection:
        # 显示示例图
        example_ge_key = 'CRISPR'
        gene_a, gene_b = 'NFKB1','RELA'
        example_l = 'Lymphoid'
        
        st.markdown("---")
        st.markdown("#### Example {}-PCorr signal of {}-{} ({})".format(example_l,gene_a,gene_b,example_ge_key))

        highlight_term = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="")

        fig = PPI_plotly(gene_a, gene_b, example_ge_key, lineage_filter=example_l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
        
    else:
        # 选中时显示对应图
        selected_row = selected_rows[0]
        ppi_str = selected_row.get("PPI", "")
        ge_key = selected_row.get("Data source", "")

        gene_a, gene_b = ppi_str.split("-")
        
        if ge_key == 'CRISPR; RNAi':
            st.markdown("---")
            st.markdown("#### {}-PCorr signal of {}-{} (CRISPR)".format(l,gene_a,gene_b))
            
            highlight_c = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="",
                                        key='highlight_c')
            
            fig = PPI_plotly(gene_a, gene_b, 'CRISPR', lineage_filter=l, highlight_by=highlight_c.strip() if highlight_c else None)
            st.plotly_chart(fig, use_container_width=False)
            
            
            st.markdown("---")
            st.markdown("#### {}-PCorr signal of {}-{} (RNAi)".format(l,gene_a,gene_b))
            
            highlight_r = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="",
                                        key='highlight_r')
            
            fig = PPI_plotly(gene_a, gene_b, 'RNAi', lineage_filter=l, highlight_by=highlight_r.strip() if highlight_r else None)
            st.plotly_chart(fig, use_container_width=False)
        else:
            st.markdown("---")
            st.markdown("#### {}-PCorr signal of {}-{} ({})".format(l,gene_a,gene_b,ge_key))
            
            highlight_term = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="")
            
            fig = PPI_plotly(gene_a, gene_b, ge_key, lineage_filter=l, highlight_by=highlight_term.strip() if highlight_term else None)
            st.plotly_chart(fig, use_container_width=False)

elif interaction_type == 'DGI':
    if not has_selection:
        
        # ---- 显示默认化合物结构 ----
        st.markdown("---")
        st.markdown("#### Example PLX-4720 structure")
        
        example_smiles = "CCCS(=O)(=O)Nc1ccc(F)c(C(=O)c2c[nH]c3ncc(Cl)cc23)c1F"
        draw_molecule(example_smiles)

        # 显示示例图
        example_dr_key = 'CTRP'
        example_ge_key = 'CRISPR'
        drugname = 'PLX-4720'
        drug = name_id[drugname]
        gene = 'BRAF'
        example_l = 'Skin'
        
        st.markdown("---")
        st.markdown("#### Example {}-PCorr signal of {}-{}".format(example_l,drugname,gene))

        highlight_term = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="")

        fig = DGI_plotly(drug, gene, example_dr_key, example_ge_key, example_l, highlight_by=highlight_term.strip() if highlight_term else None)
        st.plotly_chart(fig, use_container_width=False)
        
    else:
        # 选中时显示对应图
        selected_row = selected_rows[0]
        dgi_str = selected_row.get("DGI", "")
        dr_key = selected_row.get("Data source (Drug)", "")
        ge_key = selected_row.get("Data source (Gene)", "")
        
        drug, gene = dgi_str.split("-")
        drugname = sorted(list(id_name[drug]))[0]
        st.markdown("---")
        st.markdown("#### {}-PCorr signal of {}-{}".format(l,drugname,gene))
        
        if ge_key == 'CRISPR; RNAi':
            highlight_c = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="",
                                        key='highlight_c')
            
            fig = DGI_plotly(drug, gene, dr_key, 'CRISPR', l, highlight_by=highlight_c.strip() if highlight_c else None)
            st.plotly_chart(fig, use_container_width=False)
            
            highlight_r = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="",
                                        key='highlight_r')
            
            fig = DGI_plotly(drug, gene, dr_key, 'RNAi', l, highlight_by=highlight_r.strip() if highlight_r else None)
            st.plotly_chart(fig, use_container_width=False)
        else:
            highlight_term = st.text_input("🔍 Search for a cancer cell line (CCL) to highlight", value="")
            
            fig = DGI_plotly(drug, gene, dr_key, ge_key, l, highlight_by=highlight_term.strip() if highlight_term else None)
            st.plotly_chart(fig, use_container_width=False)


















