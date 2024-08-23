import os
import sys
import random
import pandas as pd
import streamlit as st

from utils import create_umap, image_formatter, featurize_morgan, clean_df
from plots import plot_umap

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro, library_filenames
from info import model_urls as model_urls_list


st.set_page_config(layout="wide", page_title='AI/ML DD Workshop', page_icon=':microbe:', initial_sidebar_state='collapsed')

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

@st.cache_data
def read_library(library_filename):
    print(os.path.join(data_dir, library_filename))
    return pd.read_csv(os.path.join(data_dir, library_filename), names=["SMILES"])

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
def process_umap(df):
    df["fp"] = featurize_morgan(df["SMILES"])
    df = clean_df(df)
    
    transformed_data = create_umap(df["fp"].tolist())
    df['UMAP1'] = transformed_data.T[0]
    df['UMAP2'] = transformed_data.T[1]

    df['image'] = [image_formatter(s) for s in df['SMILES']]
    df["source_filenames"] = df.index #df["ID"]
    return df

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
if cols[0].button("Upload"):
    with st.spinner():
        filenames, df_list = process_csv_files(file_list)
        st.session_state["filenames"] = filenames
        st.session_state["df_list"] = df_list

        for i, file in enumerate(filenames):
            cols[0].write(file + " length:" + str(df_list[i].shape[0]))
    

if 'chem_space_button' not in st.session_state:
    st.session_state['chem_space_button'] = True
def toggle_chem_space():
    st.session_state['chem_space_button'] = not st.session_state['chem_space_button']
if st.button('View Chemical Space', on_click=toggle_chem_space):
    # Section 2: 
    if st.session_state['chem_space_button']:
        st.divider()
        st.header("Chemical Space Plots")
        cols2 = st.columns(2)
    
        if 'umap_points' not in st.session_state:
            umap_list = []
            df_list = st.session_state["df_list"]
            for df in df_list:
                umap_list = process_umap(df)
        
            fig = plot_umap(umap_list)
            cols2[0].write(fig)
            st.session_state['umap_points'] = umap_list
            st.session_state['plot_avail'] = True
            
        elif st.session_state['chem_space_button'] and 'plot_avail' in st.session_state:
            fig = plot_umap(st.session_state['umap_points'])
            cols2[0].write(fig)
        

