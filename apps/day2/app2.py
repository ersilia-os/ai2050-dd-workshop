import os
import sys
import random
import pandas as pd
import numpy as np
import streamlit as st

from ersilia_client import ErsiliaClient

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import model_urls as model_urls_list

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

@st.cache_data
def sample_model_urls():
    random.seed(st.session_state.random_seed)
    model_urls = {
        "eos4tw0": random.choice(model_urls_list["eos4tw0"]),
        "eos4u6p": random.choice(model_urls_list["eos4u6p"]),
        "eos3804":random.choice(model_urls_list["eos3804"])
    }
    return model_urls

model_urls = sample_model_urls()

@st.cache_resource
def get_client(model_id):
    return ErsiliaClient(model_urls[model_id])
clients = {model_id: get_client(model_id) for model_id in model_urls.keys()}

@st.cache_data(show_spinner=False)
def run_predictive_models(model_ids, smiles_list):
    results = {}
    for  model_id in model_ids:
        client = clients[model_id]
        result = client.run(smiles_list)
        results[model_id] = result
    df = pd.DataFrame({"SMILES": smiles_list})
    for model_id, result in results.items():
        columns = list(result.columns)
        columns = [column for column in columns if column != "input"]
        df = pd.concat([df, result[columns]], axis=1)
    return df

smiles_list=["CCCOCCCN", "CCCCCCC"]

dp = run_predictive_models(["eos3804"], smiles_list)

st.write(dp)