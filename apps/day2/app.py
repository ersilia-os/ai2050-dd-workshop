import os
import sys
import random
import pandas as pd
import streamlit as st

from ersilia_client import ErsiliaClient

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro, exp, q1, q2, q3, library_filenames
from info import model_urls as model_urls_list

from utils import load_acinetobacter_training_data, binarize_acinetobacter_data, lolp_reducer

from plots import plot_act_inact, plot_lolp

st.set_page_config(layout="wide", page_title='AI/ML DD Workshop', page_icon=':microbe:', initial_sidebar_state='collapsed')

if "random_seed" not in st.session_state:
    st.session_state.random_seed = random.randint(0, 10000)

# ABOUT
st.sidebar.title("About")
for i in range(4):
    st.sidebar.write(about[i])
    
# MAIN
st.title(":microbe: AI2050 - AI/ML for Drug Discovery Workshop :pill:")
st.info(intro)

# Section 1: 

st.header("Step 1: Understand your data")
col1, col2 = st.columns([2, 1])
col1.write(exp)

@st.cache_data(show_spinner=False)
def load_training_data(datapath):
    return load_acinetobacter_training_data(datapath)

datapath = os.path.join(data_dir, "eos3804.csv")
df = load_training_data(datapath)
col1.write(df)

q1_header = "Let's first think about..."
col2.write(q1_header)
q_comb = '  \n'.join(q1)
col2.info(q_comb, icon=":material/quiz:")
col2.success("Read more about this dataset in its source publication [Liu et al, 2023](https://www.nature.com/articles/s41589-023-01349-8)", icon=":material/help_clinic:")

if 'step1_button' not in st.session_state:
    st.session_state['step1_button'] = False
def toggle_mock():
    st.session_state['step1_button'] = not st.session_state['step1_button']
if st.button('Move to next section', on_click=toggle_mock):
    pass
if st.session_state['step1_button']:
# Section 2: Build a classifier
    st.divider()
    st.header("Step 2: Select a cut-off to train a classifier")

    # display metrics and slider for activity cut-off
    cols = st.columns(5)
    cols[0].metric("Mean growth", round(df["Mean"].mean(), 3))
    cols[1].metric("Standard deviation", round(df["Mean"].std(), 3))
    activity_cutoff = cols[2].slider("Activity cutoff", min_value=0.1, max_value=2., value=1., step=0.001, format="%.3f")
    dt = binarize_acinetobacter_data(df, cutoff=activity_cutoff)
    cols[3].metric("Number of actives", sum(dt["Binary"]))
    cols[4].metric("Number of inactives", dt.shape[0] - sum(dt["Binary"]))

    # display data and graph according to activity cut-off
    cols = st.columns([0.45,0.275,0.275])
    dt_ = dt[["SMILES","Mean", "Binary"]]
    cols[0].write(dt_)
    dt_["Molecule index"] = dt_.index
    fig = plot_act_inact(dt_)
    cols[1].altair_chart(fig, use_container_width=True)
    q2_comb = '  \n'.join(q2)
    cols[2].info(q2_comb,icon=":material/quiz:")

    if 'cutoff' not in st.session_state:
        st.session_state['cutoff'] = None
    def set_cutoff():
        st.session_state['cutoff'] = activity_cutoff
    if st.button("Cut-off selected", on_click=set_cutoff):
        st.write(f"Cut-off value saved: {st.session_state['cutoff']}")
    if st.session_state["cutoff"]!= None:
    # Section 3: Featurisation of molecules
        st.divider()
        st.header("Step 3: Featurise molecules")

        @st.cache_data
        def sample_model_urls():
            random.seed(st.session_state.random_seed)
            model_urls = {
                "eos4tw0": random.choice(model_urls_list["eos4tw0"]),
                "eos4u6p": random.choice(model_urls_list["eos4u6p"])
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
        
        smiles_list = df["SMILES"]
        cols = st.columns([1,1,1])
        cols[0].text("1D Descriptor: Morgan Fingerpringts")
        if 'morgan_active' not in st.session_state:
            st.session_state.dp = None
            st.session_state['morgan_active'] = False
        def toggle_model_predictions_state():
            st.session_state['morgan_active'] = not st.session_state['morgan_active']
        if cols[0].button(":rocket: Calculate Morgan descriptors", on_click=toggle_model_predictions_state):
            if not st.session_state["morgan_active"]:
                pass
            else:
                with st.spinner("Running models..."):
                    #morgan = run_predictive_models(["eos4wt0"], smiles_list)
                    morgan = pd.read_csv(os.path.join(data_dir, "eos4wt0_preds.csv"))
                    cols[0].write(morgan)
                    X = morgan.iloc[:, 2:]
                    y = [1 if x <= st.session_state["cutoff"] else 0 for x in df["Mean"]]
                    st.session_state["lolp_morgan"] = lolp_reducer(X, y)
                    X_lolp = st.session_state["lolp_morgan"]["X"]
                    @st.cache_data(show_spinner=False)
                    def do_plot_lolp(X, y):
                        return plot_lolp(X, y)
                    fig1 = do_plot_lolp(X_lolp, y)
                    cols[0].altair_chart(fig1, use_container_width=True)


        cols[1].text("3D Descriptor: Chemical Checker")
        if 'cc_active' not in st.session_state:
            st.session_state.dp = None
            st.session_state['cc_active'] = False
        def toggle_model_predictions_state():
            st.session_state['cc_active'] = not st.session_state['cc_active']
        if cols[1].button(":rocket: Calculate Chemical Checker signatures", on_click=toggle_model_predictions_state):
            if not st.session_state["cc_active"]:
                pass
            else:
                with st.spinner("Running models..."):
                    #dp = run_predictive_models(["eos4u6p"], smiles_list)
                    cc  = pd.read_csv(os.path.join(data_dir, "eos4u6p_preds.csv"))
                    cols[1].write(cc)
        
        q3_comb = '  \n'.join(q3)
        cols[2].info(q3_comb,icon=":material/quiz:")