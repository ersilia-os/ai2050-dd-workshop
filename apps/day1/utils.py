import numpy as np
import pandas as pd
import os
import copy

import streamlit as st
from io import BytesIO
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw, rdFingerprintGenerator, Descriptors, QED, Crippen
from skimage import img_as_ubyte
import base64
import umap
from sklearn.decomposition import PCA

def clean_df(df):
    """
    df: a Pandas dataframe with a column "SMILES"
    """
    
    df["molecule_index"] = df.index
    df.drop_duplicates(subset ="SMILES", keep = 'first', inplace = True)
    df = df.dropna()
    df.reset_index(level=None, drop=True, inplace=True, col_level=0, col_fill='')

    return df

@st.cache_data
def featurize_morgan(smiles_list):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
    X = []
    for s in smiles_list:
        try:
            mol = Chem.MolFromSmiles(s)
            fp = mfpgen.GetFingerprint(mol)
            X.append(fp)
        except:
            X.append(None)
    return X

@st.cache_data
def process_csv_files(filename_list):
    filenames, df_list = [], []
    for i, file in enumerate(filename_list):
        try:
            filenames.append(file.name.split(".csv")[0])
        except:
            filenames.append(file[5:].split(".csv")[0])
        try:
            df = pd.read_csv(file, sep=None)
        except:
            df = pd.read_csv(file, sep=",")
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

def combine_dfs(df_list):
    combined_df = pd.DataFrame(columns=["SMILES", "fp", "file_name", "molecule_index"])
    for df in df_list:
        combined_df = pd.concat([combined_df, df[["SMILES", "fp", "file_name", "molecule_index"]]], axis=0)
    return combined_df

@st.cache_data
def create_umap(fp_list):
    umap_transformer = umap.UMAP(a=0.001, b=1.5, min_dist=0.1)
    return umap_transformer.fit_transform(fp_list)

@st.cache_data
def create_pca(fp_list):
    pca_transformer = PCA(n_components = 2)
    return pca_transformer.fit_transform(fp_list)

@st.cache_data
def image_formatter(smiles):
    mol = Chem.MolFromSmiles(smiles)
    image = Draw.MolToImage(mol, size=(150,150))
    with BytesIO() as buffer:
        image.save(buffer, 'png')
        data = base64.encodebytes(buffer.getvalue()).decode('utf-8')
    return f"data:image/png;base64,{data}"

@st.cache_data
def calc_mol_props(df):
    mols = [Chem.MolFromSmiles(smi) for smi in df["SMILES"]]
    df['mol weight'] = [Descriptors.ExactMolWt(m) for m in mols]
    df['logp'] = [Crippen.MolLogP(m) for m in mols]
    df['qed'] = [QED.default(m) for m in mols]
    return df
