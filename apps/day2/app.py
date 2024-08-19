import os
import sys
import random
import pandas as pd
import streamlit as st

from ersilia_client import ErsiliaClient

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "..", "..", "data", "day1"))

from info import about, intro, library_filenames
from info import model_urls as model_urls_list


st.set_page_config(layout="wide", page_title='AI/ML DD Workshop', page_icon=':microbe:', initial_sidebar_state='collapsed')

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

@st.cache_data
def sample_model_urls():
    random.seed(st.session_state.random_seed)
    model_urls = {
        "eos9ei3": random.choice(model_urls_list["eos9ei3"]),
        "eos43at": random.choice(model_urls_list["eos43at"])
    }
    return model_urls

model_urls = sample_model_urls()

@st.cache_resource
def get_client(model_id):
    return ErsiliaClient(model_urls[model_id])

@st.cache_data
def read_library(library_filename):
    return list(pd.read_csv(os.path.join(data_dir, library_filename))["smiles"])

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

# ABOUT
st.sidebar.title("About")
for i in range(4):
    st.sidebar.write(about[i])
    
# MAIN
st.title(":microbe: AI2050 - AI/ML for Drug Discovery Workshop :pill:")
st.markdown(intro, unsafe_allow_html=True)

# Section 1: 
st.header("Section 1")
cols = st.columns(2)
smiles = read_library(library_filenames["Example library"])
cols[0].write(smiles)
if cols[0].button("Run predictions"):
    df = run_predictive_models(["eos9ei3", "eos43at"], smiles)
    st.session_state["preds_ready"] = df
if "preds_ready" in st.session_state:
    cols[1].write(st.session_state["preds_ready"])

if 'mock_button' not in st.session_state:
    st.session_state['mock_button'] = False
def toggle_mock():
    st.session_state['mock_button'] = not st.session_state['mock_button']
if st.button('Move to next example', on_click=toggle_mock):
        pass

# Section 2: 
if st.session_state['mock_button']:
    st.divider()
    st.header("Section 2")
