import os
import numpy as np
import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from lol import LOL


def load_acinetobacter_training_data(datapath):
    df = pd.read_csv(datapath)
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    df.reset_index(drop=True, inplace=True)
    df["Mean"] = pd.to_numeric(df["Mean"])
    return df

def binarize_acinetobacter_data(data, cutoff):
    y = [1 if i <= cutoff else 0 for i in data["Mean"]]
    columns = list(data.columns)[:5]
    data = data[columns]
    data["Binary"] = y
    return data

def lolp_reducer(X, y):
    reducer = LOL(100)
    X_ = reducer.fit_transform(X, y)
    results = {
        "reducer": reducer,
        "X": X_,
        "y": y
    }
    return results

def train_acinetobacter_ml_model(X,y):
    X = np.array(X)
    y = np.array(y)
    model = RandomForestClassifier(n_jobs=-1)
    aurocs = []
    cv_data = []
    st.toast("Setting up ML model")
    for i in range(5):
        st.toast("CV iteration {0}".format(i+1))
        print("CV iteration", i)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
        model.fit(X_train, y_train)
        y_pred = model.predict_proba(X_test)[:,1]
        aurocs += [roc_auc_score(y_test, y_pred)]
        cv_data += [(y_test, y_pred)]
    st.toast("Fitting the final model")
    print("Fitting final model")
    model.fit(X, y)
    results = {
        "model": model,
        "aurocs": aurocs,
        "cv_data": cv_data,
        "X": X,
        "y": y
    }
    return results

def predict_acinetobacter_ml_model(X, model):
    y_pred = model.predict_proba(X)[:,1]
    return y_pred

def filter_valid_smiles(smiles_list):
    smiles_list = [smiles for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]
    return [smiles for smiles in smiles_list if smiles != "" and smiles is not None]

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=(200, 200))