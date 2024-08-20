import os
import sys
import random
import pandas as pd
import streamlit as st

from ersilia_client import ErsiliaClient

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro
from info import model_urls as model_urls_list
from info import seed_molecules
from utils import get_model_info_from_github, random_molecules
from plots import markdown_card, plot_single_molecule

st.set_page_config(layout="wide", page_title='AI/ML DD Workshop DAY 3', page_icon=':microbe:', initial_sidebar_state='collapsed')

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

@st.cache_data
def sample_model_urls():
    random.seed(st.session_state.random_seed)
    model_urls = {}
    for k, v in model_urls_list.items():
        if v is None:
            model_urls[k] = None
        else:
            model_urls[k] = random.choice(v)
    return model_urls

model_urls = sample_model_urls()

@st.cache_resource
def get_client(model_id):
    url = model_urls[model_id]
    if url is None:
        return None
    return ErsiliaClient(url)

clients = {model_id: get_client(model_id) for model_id in model_urls.keys()}

@st.cache_data
def get_models_info():
    models_info = {}
    for k, _ in model_urls_list.items():
        models_info[k] = get_model_info_from_github(k)
    return models_info

models_info = get_models_info()

# ABOUT
st.sidebar.title("About")
for i in range(4):
    st.sidebar.write(about[i])
    
# MAIN
st.title(":microbe: AI2050 - AI/ML for Drug Discovery Workshop :pill:")
st.info(intro)

# Step 1
st.header("Step 1: Sample the chemical space around a seed molecule")
cols = st.columns(3)
long_seed_molecules = sorted(["{0}: {1}".format(k, v) for k, v in seed_molecules.items()])
sel_smiles = cols[0].radio(label="Seed molecules", options=long_seed_molecules).split(": ")[1]
plot_single_molecule(cols[0], sel_smiles)

long_model_ids = sorted(["{0}: {1}".format(k, v["Title"]) for k, v in models_info.items()])
sel_model = cols[0].radio(label="Ersilia Model Hub identifiers", options=long_model_ids).split(":")[0]

markdown_card(cols[0], sel_model, models_info)

button_sample_chemspace = cols[0].button("Sample molecules!")

if button_sample_chemspace:
    client = clients[sel_model]
    if client is None:
        sampled_smiles = random_molecules(100)
    df = pd.DataFrame({"SMILES": sampled_smiles})
    cols[1].dataframe(df)