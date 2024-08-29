import os
import sys
import random
import pandas as pd
import streamlit as st

from utils import *
from plots import plot_umap, plot_pca

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro, questions

st.set_page_config(layout="wide", page_title='AI/ML DD Workshop', page_icon=':microbe:', initial_sidebar_state='collapsed')

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

def describe_mols(df_list, filenames):
    prepared_dfs = []
    for i, df in enumerate(df_list):
        tmp_df = process_smiles(df)
        tmp_df["file_name"] = filenames[i]
        prepared_dfs.append(tmp_df)
    return prepared_dfs    

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
st.info(intro)

# Section 1: 
st.header("Upload Data")
cols = st.columns(2)

file_list = cols[0].file_uploader("Upload chemical datasets as CSV files. Ensure they have a 'SMILES' column", accept_multiple_files=True)
if cols[0].button("Check Files and Featurize"):
    with st.spinner():
        st.toast("Loading files")
        filenames, df_list = process_csv_files(file_list)
        st.toast("Describing molecules")
        prepared_dfs = describe_mols(df_list, filenames)

        st.session_state["filenames"] = filenames
        st.session_state["raw_dfs"] = df_list
        st.session_state["prepared_dfs"] = prepared_dfs

if "filenames" in st.session_state:
    filenames = st.session_state["filenames"]
    filesizes_raw = [df.shape[0] for df in st.session_state["raw_dfs"]]
    filesizes_processed = [df.shape[0] for df in st.session_state["prepared_dfs"]]
    files_summary_df = pd.DataFrame(zip(filenames, filesizes_raw, filesizes_processed), columns=["File Name", "Total SMILES", "Featurized SMILES"])
    
    with cols[0]:
        st.subheader("File Summary")
        st.write(files_summary_df)

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

            with st.spinner():
                plot_dfs = st.session_state["prepared_dfs"]
    
                st.toast("Preparing UMAP plot")
                umap_df = process_umap(plot_dfs)
                fig_umap = plot_umap(umap_df)
    
                st.toast("Preparing PCA plot")
                pca_df = process_pca(plot_dfs)
                fig_pca = plot_pca(pca_df)
                
                cols2 = st.columns([3,3,2])
                with cols2[0]:
                    st.subheader("UMAP")
                    st.altair_chart(fig_umap)
                with cols2[1]:
                    st.subheader("PCA")
                    st.altair_chart(fig_pca)
                with cols2[2]:
                    questions_comb = '  \n'.join(questions)
                    st.info(questions_comb, icon=":material/quiz:")

