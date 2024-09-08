import os
import sys
import random
import pandas as pd
import streamlit as st

from utils import *
from plots import plot_umap, plot_pca, plot_mol_weights, plot_logp, plot_qed, plot_legend

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro, questions, library_checkbox_names, library_filenames

st.set_page_config(layout="wide", page_title='AI/ML DD Workshop', page_icon=':microbe:', initial_sidebar_state='collapsed')

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

@st.cache_data
def describe_mols(df_list, filenames):
    prepared_dfs = []
    for i, df in enumerate(df_list):
        tmp_df = process_smiles(df)
        tmp_df["file_name"] = filenames[i]
        prepared_dfs.append(tmp_df)
    return prepared_dfs    

def process_umap(combined_df):
    transformed_data = create_umap(combined_df["fp"].tolist())
    combined_df['UMAP1'] = transformed_data.T[0]
    combined_df['UMAP2'] = transformed_data.T[1]
    combined_df['image'] = [image_formatter(s) for s in combined_df['SMILES']]
    
    return combined_df

def process_pca(combined_df):
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
cols = st.columns(2)
cols[0].header("Upload Data")

file_list = cols[0].file_uploader("Upload chemical datasets as CSV files. Ensure they have a 'SMILES' column", accept_multiple_files=True)

cols[1].header("Example Libraries")
example_libraries = []
for lib in library_checkbox_names:
    example_libraries.append(cols[1].checkbox(lib))

if cols[0].button("Check Files and Featurize"):
    with st.spinner():
        st.toast("Loading files")
        for i, lib in enumerate(example_libraries):
            if lib:
                file_list.append(os.path.join("data", library_filenames[i]))
        filenames, df_list = process_csv_files(file_list)
        
        st.toast("Describing molecules")
        prepared_dfs = describe_mols(df_list, filenames)

        st.session_state["filenames"] = filenames
        st.session_state["raw_dfs"] = df_list
        st.session_state["cleaned_dfs"] = prepared_dfs

        if "prepared_dfs" in st.session_state:
            del st.session_state["prepared_dfs"]
        if "umap" in st.session_state:
            del st.session_state['umap']
            del st.session_state['pca']
            st.session_state['chem_space_button'] = False

if "cleaned_dfs" in st.session_state:
    filenames = st.session_state["filenames"]
    filesizes_raw = [df.shape[0] for df in st.session_state["raw_dfs"]]
    filesizes_processed = [df.shape[0] for df in st.session_state["cleaned_dfs"]]
    final_prepared_dfs = [df.sample(n=2000) if df.shape[0] > 2000 else df for df in st.session_state["cleaned_dfs"]]
    filesizes_final = [df.shape[0] for df in final_prepared_dfs]
    st.session_state["prepared_dfs"] = final_prepared_dfs
    files_summary_df = pd.DataFrame(zip(filenames, filesizes_raw, filesizes_processed, filesizes_final), columns=["File Name", "Total SMILES", "Featurized SMILES", "Final Number"])
    
    with cols[0]:
        st.subheader("File Summary")
    st.write(files_summary_df)

def toggle_chem_space():
    st.session_state['chem_space_button'] = not st.session_state['chem_space_button']

if 'chem_space_button' not in st.session_state:
    st.session_state['chem_space_button'] = False
if "prepared_dfs" in st.session_state:
    if st.button('View Chemical Space', on_click=toggle_chem_space) and st.session_state['chem_space_button'] == False:
        del st.session_state['umap']
        del st.session_state['pca']

if st.session_state['chem_space_button']:
    # Section 2: 
    st.divider()
    st.header("Chemical Space Plots")

    with st.spinner():
        plot_dfs = st.session_state["prepared_dfs"]
        combined_df = combine_dfs(plot_dfs)
        st.session_state["combined_df"] = combined_df

        if 'umap' in st.session_state:
            fig_umap = st.session_state['umap']
        else:
            st.toast("Preparing UMAP plot")
            umap_df = process_umap(combined_df)
            fig_umap = plot_umap(umap_df)
            st.session_state['umap'] = fig_umap
        
        if 'pca' in st.session_state:
            fig_pca = st.session_state['pca']
        else:
            st.toast("Preparing PCA plot")
            pca_df = process_pca(combined_df)
            fig_pca = plot_pca(pca_df)
            st.session_state['pca'] = fig_pca
        
        cols2 = st.columns([3,3,2])
        with cols2[0]:
            st.subheader("UMAP")
            st.altair_chart(fig_umap)
        with cols2[1]:
            st.subheader("PCA")
            st.altair_chart(fig_pca)
        with cols2[2]:
            fig_legend = plot_legend(combined_df)
            st.write("#")
            st.altair_chart(fig_legend)
            questions_comb = '  \n'.join(questions)
            st.info(questions_comb, icon=":material/quiz:")


def toggle_chem_prop():
    st.session_state['chem_prop_button'] = not st.session_state['chem_prop_button']
    
if 'chem_prop_button' not in st.session_state:
    st.session_state['chem_prop_button'] = False
if st.session_state['chem_space_button']:
    if st.button('Plot Chemical Properties', on_click=toggle_chem_prop) and st.session_state['chem_prop_button']:
        # Section 3
        st.divider()
        st.header("Distributions of Chemical Properties")
        cols3 = st.columns([2,2,2,1])
        
        combined_df = st.session_state["combined_df"]
        combined_props = calc_mol_props(combined_df)
        
        fig_mw = plot_mol_weights(combined_props)
        cols3[0].altair_chart(fig_mw)
        fig_logp = plot_logp(combined_props)
        cols3[1].altair_chart(fig_logp)
        fig_qed = plot_qed(combined_props)
        cols3[2].altair_chart(fig_qed)
        fig_legend = plot_legend(combined_props)
        cols3[3].altair_chart(fig_legend)
        