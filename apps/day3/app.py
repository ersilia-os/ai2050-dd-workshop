# Imports
import os
import sys
import streamlit as st

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)

from components import general_setup, header_and_sidebar, select_seed_molecule, select_model_sampler, step_1_blurbs
from components import sample_molecules_from_selected_model, sampled_molecules_explorer, wishlist, prediction_model_selector, step_2_blurbs
from components import calculate_properties, explore_predictions, step_3_blurbs

if 'sampled_molecules' not in st.session_state:
    st.sampled_molecules = []

# Header and sidebar
header_and_sidebar()

# General setup
results = general_setup()
model_urls = results["model_urls"]
activity_models_urls = results["activity_models_urls"]
adme_models_urls = results["adme_models_urls"]
clients = results["clients"]
activity_clients = results["activity_clients"]
adme_clients = results["adme_clients"]
models_info = results["models_info"]
activity_models_info = results["activity_models_info"]
adme_models_info = results["adme_models_info"]

# Select seed molecule
cols = st.columns(4)
sel_smiles, sel_inchikey = select_seed_molecule(cols[0])

# Select model sampler
sel_model = select_model_sampler(cols[1], cols[2], models_info)

# Step 1 blurbs
sample_button = step_1_blurbs(cols[3])
if sample_button:
    # Sampled molecules
    df = sample_molecules_from_selected_model(sel_model, sel_inchikey, sel_smiles, clients)

# Explore sampled molecules and get the wishlist
cols = st.columns([3, 1, 3])

sampled_molecules_explorer(cols[0], cols[1], df, sel_model, models_info)

dl = wishlist(cols[2], sel_smiles, sel_inchikey)

cols = st.columns(4)
results = prediction_model_selector(cols[0], cols[1], cols[2], activity_models_urls, activity_models_info, adme_models_urls, adme_models_info)
sel_activity_model = results["sel_activity_model"]
sel_adme_model = results["sel_adme_model"]
sel_adme_prop_1 = results["sel_adme_prop_1"]
sel_adme_prop_2 = results["sel_adme_prop_2"]
sel_adme_prop_3 = results["sel_adme_prop_3"]
sel_adme_prop_4 = results["sel_adme_prop_4"]
sel_adme_prop_5 = results["sel_adme_prop_5"]

# Step 2 blurbs
step_2_blurbs(cols[3])

# Calculate properties
dl, activity_column = calculate_properties(dl, sel_inchikey, activity_clients, sel_activity_model, adme_clients, sel_adme_model, sel_adme_prop_1, sel_adme_prop_2, sel_adme_prop_3, sel_adme_prop_4, sel_adme_prop_5)

# Predictions results
cols = st.columns([1/2, 1/4, 1/4, 1/4])
explore_predictions(cols[0], cols[1], cols[2], dl, activity_column)

# Step 3 blurbs
step_3_blurbs(cols[3])


