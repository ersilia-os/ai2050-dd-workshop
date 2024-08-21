import os
import sys
import random
import pandas as pd
import streamlit as st

from ersilia_client import ErsiliaClient

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro, step_1_explanation, step_2_explanation, step_3_explanation, adme_warning_message
from info import model_urls as model_urls_list
from info import activity_models_urls as activity_models_urls_list
from info import adme_models_urls as adme_models_urls_list
from info import seed_molecules
from info import adme_col_new2old
from utils import get_model_info_from_github, random_molecules, basic_molecules_dataframe
from plots import markdown_card, plot_single_molecule, basic_mols2grid_plot

st.set_page_config(layout="wide", page_title='AI/ML DD Workshop Day 3', page_icon=':microbe:', initial_sidebar_state='collapsed')

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

@st.cache_data
def sample_activity_models_urls():
    random.seed(st.session_state.random_seed)
    activity_models_urls = {}
    for k, v in activity_models_urls_list.items():
        if v is None:
            activity_models_urls[k] = None
        else:
            activity_models_urls[k] = random.choice(v)
    return activity_models_urls

@st.cache_data
def sample_adme_model_urls():
    random.seed(st.session_state.random_seed)
    adme_models_urls = {}
    for k, v in adme_models_urls_list.items():
        if v is None:
            adme_models_urls[k] = None
        else:
            adme_models_urls[k] = random.choice(v)
    return adme_models_urls

model_urls = sample_model_urls()
activity_models_urls = sample_activity_models_urls()
adme_models_urls = sample_adme_model_urls()

@st.cache_resource
def get_client(model_id):
    url = model_urls[model_id]
    if url is None:
        return None
    return ErsiliaClient(url)

@st.cache_resource
def get_activity_client(model_id):
    url = activity_models_urls[model_id]
    if url is None:
        return None
    return ErsiliaClient(url)

@st.cache_resource
def get_adme_client(model_id):
    url = adme_models_urls[model_id]
    if url is None:
        return None
    return ErsiliaClient

clients = {model_id: get_client(model_id) for model_id in model_urls.keys()}
activity_clients = {model_id: get_activity_client(model_id) for model_id in activity_models_urls.keys()}
adme_clients = {model_id: get_adme_client(model_id) for model_id in adme_models_urls.keys()}

@st.cache_data
def get_models_info():
    models_info = {}
    for k, _ in model_urls_list.items():
        models_info[k] = get_model_info_from_github(k)
    return models_info

@st.cache_data
def get_activity_models_info():
    models_info = {}
    for k, _ in activity_models_urls_list.items():
        models_info[k] = get_model_info_from_github(k)
    return models_info

@st.cache_data
def get_adme_models_info():
    models_info = {}
    for k, _ in adme_models_urls_list.items():
        models_info[k] = get_model_info_from_github(k)
    return models_info

models_info = get_models_info()
activity_models_info = get_activity_models_info()
adme_models_info = get_adme_models_info()

# ABOUT
st.sidebar.title("About")
for i in range(4):
    st.sidebar.write(about[i])
    
# MAIN
st.title(":microbe: AI2050 - AI/ML for Drug Discovery Workshop :pill:")
st.info(intro)

# Step 1
st.header("Sample the chemical space around a seed molecule")
cols = st.columns(4)
long_seed_molecules = sorted(["{0}: {1}".format(k, v) for k, v in seed_molecules.items()])
cols[0].subheader("Select a seed molecule")
sel_smiles = cols[0].radio(label="Seed molecules", options=long_seed_molecules).split(": ")[1]
plot_single_molecule(cols[0], sel_smiles)

long_model_ids = sorted(["{0}: {1}".format(k, v["Title"]) for k, v in models_info.items()])
cols[1].subheader("Select a model")
sel_model = cols[1].radio(label="Ersilia Model Hub identifiers", options=long_model_ids).split(":")[0]

markdown_card(cols[2], sel_model, models_info)

cols[3].info(step_1_explanation)

client = clients[sel_model]
if client is None:
    sampled_smiles = random_molecules(100)


df = basic_molecules_dataframe(sampled_smiles, sel_smiles)

basic_mols2grid_plot(df)

cols = st.columns(2)

cols[0].subheader("Molecules sampled by model {0}: {1}".format(sel_model, models_info[sel_model]["Title"]))
cols[0].dataframe(df)

cols[1].subheader("Your molecules wishlist")
smiles_library = cols[1].text_area(label="Chosen molecules (SMILES)", value=sel_smiles, height=600)

dl = basic_molecules_dataframe(smiles_library.split("\n"), sel_smiles)

st.subheader("Your molecules wishlist ({0})".format(dl.shape[0]))
basic_mols2grid_plot(dl)

st.header(":pill: Predict properties of your molecules wishlist")

# Select models
cols = st.columns(4)
sel_activity_model = cols[0].radio(label="Select an activity model", options=list(activity_models_urls.keys()), index=0)
sel_adme_model = [k for k in adme_models_urls.keys()][0]
adme_columns = list(adme_col_new2old.keys())
_adme_columns = [k for k in adme_columns]
sel_adme_prop_1 = cols[0].selectbox(label="Select ADME property 1", options=_adme_columns, index=0)
_adme_columns.remove(sel_adme_prop_1)
sel_adme_prop_2 = cols[0].selectbox(label="Select ADME property 2", options=_adme_columns, index=0)
_adme_columns.remove(sel_adme_prop_2)
sel_adme_prop_3 = cols[0].selectbox(label="Select ADME property 3", options=_adme_columns, index=0)
_adme_columns.remove(sel_adme_prop_3)
sel_adme_prop_4 = cols[0].selectbox(label="Select ADME property 4", options=_adme_columns, index=0)
_adme_columns.remove(sel_adme_prop_4)
sel_adme_prop_5 = cols[0].selectbox(label="Select ADME property 5", options=_adme_columns, index=0)
_adme_columns.remove(sel_adme_prop_5)
cols[0].warning(adme_warning_message)

markdown_card(cols[1], sel_activity_model, activity_models_info)
markdown_card(cols[2], sel_adme_model, adme_models_info)

cols[3].info(step_2_explanation)


# Calculate activity
activity_client = activity_clients[sel_activity_model]

def random_activity(n):
    return [round(random.uniform(0, 1), 2) for _ in range(n)]

if activity_client is None:
    activity_scores = random_activity(dl.shape[0])

# Calculate ADME properties
adme_client = adme_clients[sel_adme_model]

def random_adme_percentiles(n, m):
    R = []
    for _ in range(n):
        r = []
        for _ in range(m):
            r.append(round(random.uniform(0, 1)*100, 2))
        R.append(r)
    return R

if adme_client is None:
    adme_percentiles = random_adme_percentiles(dl.shape[0], 5)

dl["Activity"] = activity_scores
da = pd.DataFrame(adme_percentiles, columns=[sel_adme_prop_1, sel_adme_prop_2, sel_adme_prop_3, sel_adme_prop_4, sel_adme_prop_5])
dl = pd.concat([dl[["InChIKey", "SMILES", "Tanimoto Coeff", "Activity"]], da], axis=1)

st.subheader(":star: Here are your molecules with some predicted properties!")
cols = st.columns(2)

cols[0].dataframe(dl)
cols[-1].info(step_3_explanation)

