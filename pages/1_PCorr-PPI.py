# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 13:51:09 2025

@author: zengt
"""

import streamlit as st
import pickle
import os
import pandas as pd

def table_showing():
    script_dir = os.path.dirname(os.path.dirname(__file__))
    
    # 读取数据
    data_path = os.path.join(script_dir, "data/de novo PPIs.pkl")
    df = pickle.load(open(data_path, 'rb'))
    
    column_to_filter = st.selectbox("选择要筛选的列", df.columns)
    
    # 检查所选列是否为连续性变量（数值类型）
    if pd.api.types.is_numeric_dtype(df[column_to_filter]):
        # 获取该列的最小值和最大值
        min_value = float(df[column_to_filter].min())
        max_value = float(df[column_to_filter].max())
    
        # 创建区间滑块
        selected_range = st.slider(
            f"选择 {column_to_filter} 的范围",
            min_value=min_value,
            max_value=max_value,
            value=(min_value, max_value)
        )
    
        # 筛选 DataFrame
        filtered_df = df[(df[column_to_filter] >= selected_range[0]) & (df[column_to_filter] <= selected_range[1])]
    else:
        # 如果不是数值类型，提供一个输入框让用户输入筛选值
        selected_value = st.text_input(f"输入要筛选的 {column_to_filter} 值")
        if selected_value:
            # 筛选 DataFrame
            filtered_df = df[df[column_to_filter].astype(str).str.contains(selected_value, case=False, na=False)]
        else:
            filtered_df = df
    
    # 显示筛选后的 DataFrame
    st.write("筛选后的数据：")
    st.write(filtered_df)
        
    # columns = table.columns.tolist()
    
    
    # # 创建选择框让用户选择两列
    # col1 = st.selectbox("Select Gene A", columns)
    
    # # 创建筛选条件
    # st.write("### Filter Data")
    # filter_col1 = st.selectbox(f"Select {col1} value", table[col1].unique())
    
    # # 筛选 DataFrame
    # filtered_table = table[table[col1] == filter_col1]
    
    # # 显示筛选后的 DataFrame
    # st.write("### Filtered DataFrame")
    # st.write(filtered_table)

st.set_page_config(page_title="Functional protein-protein interactions (PPIs)")


table_showing()


