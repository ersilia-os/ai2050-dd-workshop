import os
import sys
import streamlit as st
import random
import pandas as pd
from ersilia_client import ErsiliaClient

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro
from info import model_urls as model_urls_list
from info import activity_models_urls as activity_models_urls_list
from info import adme_models_urls as adme_models_urls_list
from info import seed_molecules
from info import adme_col_new2old, adme_warning_message
from info import step_1_explanation, step_1_questions
from info import step_2_explanation, step_2_questions
from info import step_3_explanation, step_3_questions

from utils import get_model_info_from_github
from utils import smiles_to_inchikey
from utils import random_molecules, cached_molecules
from utils import basic_molecules_dataframe
from utils import filter_smiles_library
from utils import cached_predictions, random_activity, random_adme_percentiles

from plots import plot_single_molecule
from plots import markdown_card
from plots import basic_mols2grid_plot
from plots import basic_molecule_card
from plots import molecule_card
from plots import histogram


def header_and_sidebar():
    st.set_page_config(layout="wide", page_title='AI/ML DD Workshop Day 3', page_icon=':microbe:', initial_sidebar_state='collapsed')
    if "random_seed" not in st.session_state:
        st.session_state.random_seed = random.randint(0, 10000)
    st.sidebar.title("About")
    for i in range(4):
        st.sidebar.write(about[i])
    st.title(":microbe: AI2050 - AI/ML for Drug Discovery Workshop :pill:")
    st.info(intro)


def general_setup():
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

    result = {
        "model_urls": model_urls,
        "activity_models_urls": activity_models_urls,
        "adme_models_urls": adme_models_urls,
        "clients": clients,
        "activity_clients": activity_clients,
        "adme_clients": adme_clients,
        "models_info": models_info,
        "activity_models_info": activity_models_info,
        "adme_models_info": adme_models_info
    }

    return result


def select_seed_molecule(col):
    long_seed_molecules = sorted(["`{0}`: {1}".format(k, v) for k, v in seed_molecules.items()])
    col.subheader("Select a seed molecule")
    sel_smiles = col.radio(label="Seed molecules", options=long_seed_molecules).split(": ")[1]
    sel_inchikey = smiles_to_inchikey(sel_smiles)
    plot_single_molecule(col, sel_smiles)
    return sel_smiles, sel_inchikey


def select_model_sampler(col_0, col_1, models_info):
    long_model_ids = sorted(["`{0}`: {1}".format(k, v["Title"]) for k, v in models_info.items()])
    col_0.subheader("Select a model")
    sel_model = col_0.radio(label="Ersilia Model Hub identifiers", options=long_model_ids).split("`:")[0].split("`")[1]
    markdown_card(col_1, sel_model, models_info)
    return sel_model


def step_1_blurbs(col):
    col.success(step_1_explanation)
    col.info(step_1_questions)
    button = col.button("Sample molecules! :crystal_ball:", key="sample_molecules")
    return button


def sample_molecules_from_selected_model(sel_model, sel_inchikey, sel_smiles, clients):
    client = clients[sel_model]
    if client is None:
        sampled_smiles = cached_molecules(inchikey=sel_inchikey, model_id=sel_model, n=100)
        if sampled_smiles is None:
            sampled_smiles = random_molecules(100)
    df = basic_molecules_dataframe(sampled_smiles, sel_smiles)
    return df


def sampled_molecules_explorer(col_0, col_1, df, sel_model, models_info):
    col_0.subheader("Molecules sampled by model {0}: {1}".format(sel_model, models_info[sel_model]["Title"]))
    scols = col_0.columns(4)
    min_tanimoto, max_tanimoto = scols[0].slider("Tanimoto Coeff", min_value=0.0, max_value=1.0, value=(0.0, 1.0))
    df_ = df[(df["Tanimoto Coeff"] >= min_tanimoto) & (df["Tanimoto Coeff"] <= max_tanimoto)]
    min_mw, max_mw = scols[1].slider("MolWeight", min_value=0, max_value=1000, value=(0, 1000))
    df_ = df_[(df_["MolWeight"] >= min_mw) & (df_["MolWeight"] <= max_mw)]
    min_logp, max_logp = scols[2].slider("LogP", min_value=-10.0, max_value=10.0, value=(-10.0, 10.0))
    df_ = df_[(df_["LogP"] >= min_logp) & (df_["LogP"] <= max_logp)]
    min_qed, max_qed = scols[3].slider("QED", min_value=0.0, max_value=1.0, value=(0.0, 1.0))
    df_ = df_[(df_["QED"] >= min_qed) & (df_["QED"] <= max_qed)]
    col_0.dataframe(df_, height=600)
    if df.shape[0] > 0:
        col_1.subheader("")
        molecule_index_to_show = col_1.number_input("Index of the molecule to explore", min_value=0, max_value=df.shape[0]-1, value=0)
        molecule_to_show = df.iloc[molecule_index_to_show]["SMILES"]
        basic_molecule_card(col_1, molecule_to_show, df)


def wishlist(col, sel_smiles, sel_inchikey):
    col.subheader("Your molecules wishlist")
    smiles_library = col.text_area(label="Chosen molecules (SMILES)", value="", height=600)
    smiles_library = [x.strip() for x in smiles_library.split("\n")]
    smiles_library = [x for x in smiles_library if x != ""]
    n = len(smiles_library)
    smiles_library = list(set(smiles_library))
    smiles_library = filter_smiles_library(sel_inchikey, smiles_library)
    if len(smiles_library) == 0:
        st.error("You have entered {0} input molecules of which {1} are unique and present in the sampled chemical space of your seed.".format(n, len(smiles_library)))
    elif n == len(smiles_library):
        st.success("You have entered {0} input molecules of which {1} are unique and present in the sampled chemical space of your seed.".format(n, len(smiles_library)))
    else:
        st.warning("You have entered {0} input molecules of which {1} are unique and present in the sampled chemical space of your seed.".format(n, len(smiles_library)))
    dl = basic_molecules_dataframe(smiles_library, sel_smiles)
    if dl.shape[0] == 0:
        st.info("Enter at least one valid molecule! Copy-paste SMILES strings from the table on the left.")
        st.stop()
        return None
    if dl.shape[0] > 1000:
        st.warning("You have entered more than 1000 molecules ({0}). We will only show the first 1000.".format(dl.shape[0]))
        dl = dl.head(1000)
    st.subheader("Your molecules wishlist ({0})".format(dl.shape[0]))
    basic_mols2grid_plot(dl)
    return dl


def prediction_model_selector(col_0, col_1, col_2, activity_models_urls, activity_models_info, adme_models_urls, adme_models_info):
    activity_model_keys = sorted(activity_models_urls.keys())
    activity_model_labels = ["`{0}`: {1}".format(k, activity_models_info[k]["Title"]) for k in activity_model_keys]
    sel_activity_model = col_0.radio(label="Select an activity model", options=activity_model_labels, index=0).split("`:")[0].split("`")[1]
    sel_adme_model = [k for k in adme_models_urls.keys()][0]
    adme_columns = list(adme_col_new2old.keys())
    _adme_columns = [k for k in adme_columns]
    sel_adme_prop_1 = col_0.selectbox(label="Select ADME property 1", options=_adme_columns, index=0)
    _adme_columns.remove(sel_adme_prop_1)
    sel_adme_prop_2 = col_0.selectbox(label="Select ADME property 2", options=_adme_columns, index=0)
    _adme_columns.remove(sel_adme_prop_2)
    sel_adme_prop_3 = col_0.selectbox(label="Select ADME property 3", options=_adme_columns, index=0)
    _adme_columns.remove(sel_adme_prop_3)
    sel_adme_prop_4 = col_0.selectbox(label="Select ADME property 4", options=_adme_columns, index=0)
    _adme_columns.remove(sel_adme_prop_4)
    sel_adme_prop_5 = col_0.selectbox(label="Select ADME property 5", options=_adme_columns, index=0)
    _adme_columns.remove(sel_adme_prop_5)
    col_0.warning(adme_warning_message)
    markdown_card(col_1, sel_activity_model, activity_models_info)
    markdown_card(col_2, sel_adme_model, adme_models_info)
    results = {
        "sel_activity_model": sel_activity_model,
        "sel_adme_model": sel_adme_model,
        "sel_adme_prop_1": sel_adme_prop_1,
        "sel_adme_prop_2": sel_adme_prop_2,
        "sel_adme_prop_3": sel_adme_prop_3,
        "sel_adme_prop_4": sel_adme_prop_4,
        "sel_adme_prop_5": sel_adme_prop_5
    }
    return results


def step_2_blurbs(col):
    col.success(step_2_explanation)
    col.info(step_2_questions)
    button = col.button("Calculate properties! :rocket:")
    return button


def calculate_properties(dl, sel_inchikey, activity_clients, sel_activity_model, adme_clients, sel_adme_model, sel_adme_prop_1, sel_adme_prop_2, sel_adme_prop_3, sel_adme_prop_4, sel_adme_prop_5):
    activity_client = activity_clients[sel_activity_model]
    if activity_client is None:
        activity_scores = cached_predictions(sel_inchikey, sel_activity_model, dl["SMILES"].tolist())
        if activity_scores is None:
            activity_column = "Activity"
            activity_scores = random_activity(dl.shape[0])
        else:
            activity_column = list(activity_scores.columns)[-1].split("_")[-1]
            activity_scores = activity_scores[list(activity_scores.columns)[-1]].tolist()
    adme_client = adme_clients[sel_adme_model]
    if adme_client is None:
        adme_percentiles = cached_predictions(sel_inchikey, sel_adme_model, dl["SMILES"].tolist())
        if adme_percentiles is None:
            adme_percentiles = random_adme_percentiles(dl.shape[0], 5)
            da = pd.DataFrame(adme_percentiles, columns=[sel_adme_prop_1, sel_adme_prop_2, sel_adme_prop_3, sel_adme_prop_4, sel_adme_prop_5])
        else:
            cols = [k for k in list(adme_percentiles.columns) if k.lower() in ["inchikey", "smiles"] or k.endswith("approved_percentile")]
            adme_percentiles = adme_percentiles[cols]
            rename = dict((k, "_".join(k.split("_")[1:])) for k in adme_percentiles.columns if k.lower() not in ["inchikey", "smiles"])
            adme_percentiles = adme_percentiles.rename(columns=rename, inplace=False)
            rename = dict((v, k) for k, v in adme_col_new2old.items())
            adme_percentiles = adme_percentiles.rename(columns=rename, inplace=False)
            da = adme_percentiles[[sel_adme_prop_1, sel_adme_prop_2, sel_adme_prop_3, sel_adme_prop_4, sel_adme_prop_5]]
    dl[activity_column] = activity_scores
    dl = pd.concat([dl[["InChIKey", "SMILES", "Tanimoto Coeff", activity_column]], da], axis=1)
    return dl, activity_column


def explore_predictions(col_0, col_1, col_2, dl, activity_column):
    col_0.subheader(":star: Here are your molecules with some predicted properties!")
    min_activity, max_activity = col_0.slider(activity_column, min_value=0.0, max_value=1.0, value=(0.0, 1.0))
    dl_ = dl[(dl[activity_column] >= min_activity) & (dl[activity_column] <= max_activity)]
    col_0.dataframe(dl_, height=600)
    col_0.success("These are your selected molecules with their calculated properties. Feel free to explore them further!")
    molecule_index_to_show = col_1.number_input("Index of the molecule to select", min_value=0, max_value=dl.shape[0]-1, value=0)
    explore_smiles = dl.iloc[molecule_index_to_show]["SMILES"]
    molecule_card(col_1, explore_smiles, dl, activity_column)
    column = col_2.selectbox("Explore the distribution of column values", options=list(dl.columns)[2:])
    histogram(col_2, dl, column)


def step_3_blurbs(col):
    col.success(step_3_explanation)
    col.info(step_3_questions)