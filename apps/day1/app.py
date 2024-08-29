import os
import sys
import random
import pandas as pd
import streamlit as st
import copy

from utils import create_umap, create_pca, image_formatter, featurize_morgan, clean_df
from plots import plot_umap, plot_pca

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro, library_filenames

st.set_page_config(layout="wide", page_title='AI/ML DD Workshop', page_icon=':microbe:', initial_sidebar_state='collapsed')

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

def process_csv_files(filename_list):
    filenames, df_list = [], []
    for i, file in enumerate(filename_list):
        filenames.append(file.name.split(".csv")[0])
        df = pd.read_csv(file, sep=None)
        for col in df.columns:
            if "SMILES" in col.upper():
                df.rename(columns = {col : "SMILES"}, inplace=True)
        df = df.loc[:,~df.columns.duplicated()].copy()
        df_list.append(df)
    return filenames, df_list
        
@st.cache_data
def process_smiles(df):
    tmp_df = copy.copy(df)
    tmp_df["fp"] = featurize_morgan(tmp_df["SMILES"])
    tmp_df = clean_df(tmp_df[["SMILES", "fp"]])
    return tmp_df

def process_umap(df_list):
    combined_df = pd.DataFrame(columns=["SMILES", "fp", "file_name", "molecule_index"])
    for df in df_list:
        combined_df = pd.concat([combined_df, df[["SMILES", "fp", "file_name", "molecule_index"]]], axis=0)
    
    transformed_data = create_umap(combined_df["fp"].tolist())
    combined_df['UMAP1'] = transformed_data.T[0]
    combined_df['UMAP2'] = transformed_data.T[1]
    combined_df['image'] = [image_formatter(s) for s in combined_df['SMILES']]
    
    return combined_df

def process_pca(df_list):
    combined_df = pd.DataFrame(columns=["SMILES", "fp", "file_name", "molecule_index"])
    for df in df_list:
        combined_df = pd.concat([combined_df, df[["SMILES", "fp", "file_name", "molecule_index"]]], axis=0)
    
    transformed_data = create_pca(combined_df["fp"].tolist())
    combined_df['PCA1'] = transformed_data.T[0]
    combined_df['PCA2'] = transformed_data.T[1]
    combined_df['image'] = [image_formatter(s) for s in combined_df['SMILES']]
    
    return combined_df

# ABOUT
st.sidebar.title("About")
for i in range(4):
    st.sidebar.write(about[i])
    
# MAIN
st.title(":microbe: AI2050 - AI/ML for Drug Discovery Workshop :pill:")
st.markdown(intro, unsafe_allow_html=True)

# Section 1: 
st.header("Upload Data")
cols = st.columns(2)

file_list = cols[0].file_uploader("Upload chemical datasets as CSV files. Ensure they have a 'SMILES' column", accept_multiple_files=True)
if cols[0].button("Check Files"):
    with st.spinner():
        filenames, df_list = process_csv_files(file_list)
        st.session_state["filenames"] = filenames
        st.session_state["df_list"] = df_list
        st.session_state['chem_space_button'] = False

        filesizes = [df.shape[0] for df in df_list]
        files_df = pd.DataFrame(zip(filenames, filesizes), columns=["File Name", "Total SMILES"])
        with cols[0]:
            st.subheader("File Summary")
            st.write(files_df)

def toggle_chem_space():
    st.session_state['chem_space_button'] = not st.session_state['chem_space_button']

if 'chem_space_button' not in st.session_state:
    st.session_state['chem_space_button'] = False
if "filenames" in st.session_state:
    if st.button('View Chemical Space', on_click=toggle_chem_space) and st.session_state['chem_space_button']:
        # Section 2: 
        if st.session_state['chem_space_button']:
            st.divider()
            st.header("Chemical Space Plots")

            df_list = st.session_state["df_list"]
            plot_dfs = []
            for i, df in enumerate(df_list):
                tmp_df = process_smiles(df)
                tmp_df["file_name"] = st.session_state["filenames"][i]
                plot_dfs.append(tmp_df)
            
            umap_df = process_umap(plot_dfs)
            fig_umap = plot_umap(umap_df)

            pca_df = process_pca(plot_dfs)
            fig_pca = plot_pca(pca_df)
            
            cols2 = st.columns(2)
            with cols2[0]:
                st.subheader("UMAP")
                st.altair_chart(fig_umap)
            with cols2[1]:
                st.subheader("PCA")
                st.altair_chart(fig_pca)            
            


        

