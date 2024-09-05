import os
import sys
import random
import pandas as pd
import numpy as np
import requests
from sklearn import metrics
import altair as alt
import streamlit as st

from ersilia_client import ErsiliaClient

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
data_dir = os.path.abspath(os.path.join(root, "data"))

from info import about, intro, exp, q1, q2, q3, q4, q4_followup, q5, library_filenames
from info import model_urls as model_urls_list

from utils import load_acinetobacter_training_data, binarize_acinetobacter_data, lolp_reducer, train_acinetobacter_ml_model, predict_acinetobacter_ml_model 
from utils import draw_molecule, calculate_sensitivity_recall

from plots import plot_act_inact, plot_lolp, plot_roc_curve, plot_contingency_table

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
col2.success("Read more about this dataset in its source publication [Liu et al, 2023](https://drive.google.com/file/d/1_EwSEDzbjtqpLMGA8EIrxkxxJRsTrE0r/view?usp=sharing)", icon=":material/help_clinic:")

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
        y = [1 if x <= st.session_state["cutoff"] else 0 for x in df["Mean"]] # INSTANTIATE ONLY ONE TIME; DO NOT SUBSCRIBE!
    # Section 3: Featurisation of molecules
        st.divider()
        st.header("Step 3: Featurise molecules")

        @st.cache_data(show_spinner=False)
        def load_dataframe(url, filename):
            if not os.path.exists(filename):
                response = requests.get(url)
                response.raise_for_status() 
                print(response)
                with open(filename, 'wb') as f:
                    f.write(response.content)
            print(filename)
            df = pd.read_csv(filename)
            print(df.head())
            return df
    
        @st.cache_data(show_spinner=False)
        def do_plot_lolp(X, y):
            return plot_lolp(X, y)
        
        smiles_list = df["SMILES"]
        cols = st.columns([0.35,0.35,0.3])
        q3_comb = '  \n'.join(q3)
        cols[2].info(q3_comb,icon=":material/quiz:")

        if 'desc1_active' not in st.session_state:
            st.session_state['desc1_active'] = False
        def toggle_desc1_state():
            st.session_state['desc1_active'] = not st.session_state['desc1_active']
        if 'desc1_results' not in st.session_state:
            st.session_state['desc1_results'] = None
        if 'desc1_lolp' not in st.session_state:
            st.session_state['desc1_lolp'] = None
        if cols[0].button(":rocket: Calculate Morgan descriptors", on_click=toggle_desc1_state):
            if not st.session_state["desc1_active"]:
                pass
            else:
                with st.spinner("Running Ersilia model..."):
                    url = "https://ai2050-workshops.s3.eu-central-1.amazonaws.com/eos4wt0_preds.csv"
                    desc1 = load_dataframe(url, "eos4wt0_preds.csv")
                    st.session_state['desc1_results'] = desc1
                    X = desc1.iloc[:, 2:]
                    st.session_state["desc1_lolp"] = lolp_reducer(X, y)
        if st.session_state['desc1_results'] is not None:
            cols[0].write(st.session_state['desc1_results'])
            fig1 = do_plot_lolp(st.session_state["desc1_lolp"]["X"], y)
            cols[1].altair_chart(fig1, use_container_width=True)

        cols = st.columns([0.35,0.35,0.3])
        if 'desc2_active' not in st.session_state:
            st.session_state['desc2_active'] = False
        def toggle_desc2_state():
            st.session_state['desc2_active'] = not st.session_state['desc2_active']
        if 'desc2_results' not in st.session_state:
            st.session_state['desc2_results'] = None
        if 'desc2_lolp' not in st.session_state:
            st.session_state['desc2_lolp'] = None
        if cols[0].button(":rocket: Calculate Chemical Checker signatures", on_click=toggle_desc2_state):
            if not st.session_state["desc2_active"]:
                pass
            else:
                with st.spinner("Running Ersilia model..."):
                    #url = "https://ai2050-workshops.s3.eu-central-1.amazonaws.com/eos4u6p_preds.csv"
                    url = "https://ai2050-workshops.s3.eu-central-1.amazonaws.com/eos4u6p_preds_red.csv"
                    desc2 = load_dataframe(url, "eos4u6p_preds.csv")
                    st.session_state['desc2_results'] = desc2
                    X = desc2.iloc[:, 2:]
                    st.session_state['desc2_lolp'] = lolp_reducer(X, y)
        if st.session_state['desc2_results'] is not None:
            cols[0].write(st.session_state['desc2_results'])
            fig1 = do_plot_lolp(st.session_state["desc2_lolp"]["X"], y)
            cols[1].altair_chart(fig1, use_container_width=True)

            st.divider()
            st.header("Step 4: Train a classifier")
            cols = st.columns([0.3, 0.35, 0.35])
            q4_comb = '  \n'.join(q4)
            cols[0].info(q4_comb,icon=":material/quiz:")

            @st.cache_data(show_spinner=False)
            def do_plot_roc_curve(tprs_df):
                return plot_roc_curve(tprs_df)

            if 'train_desc1_model_active' not in st.session_state:
                st.session_state['train_desc1_model_active'] = False
            def toggle_train_desc1_model_state():
                st.session_state['train_desc1_model_active'] = not st.session_state['train_desc1_model_active']
            if cols[1].button('🤖 Train a classifier using Morgan FPS!', on_click=toggle_train_desc1_model_state):
                if not st.session_state["train_desc1_model_active"]:
                    pass
                else:
                    with st.spinner("Training the model..."):
                        st.session_state.model_results_desc1 = train_acinetobacter_ml_model(st.session_state['desc1_results'].iloc[:, 2:], y)
            if st.session_state["train_desc1_model_active"]:
                if "model_results_desc1" in st.session_state:
                    aurocs = st.session_state.model_results_desc1["aurocs"]
                    std_auroc = np.std(aurocs)
                    tprs = []
                    mean_fpr = np.linspace(0,1,100)
                    for i in st.session_state.model_results_desc1["cv_data"]:
                        fpr, tpr, _ = metrics.roc_curve(i[0],i[1])
                        interp_tpr = np.interp(mean_fpr, fpr, tpr)
                        interp_tpr[0] = 0.0
                        tprs.append(interp_tpr)
                    mean_tpr = np.mean(tprs, axis=0)
                    mean_tpr[-1] = 1.0
                tprs_df = pd.DataFrame({
                    'tpr_cv1': tprs[0],
                    'tpr_cv2': tprs[1],
                    'tpr_cv3': tprs[2],
                    'tpr_cv4': tprs[3],
                    'tpr_cv5': tprs[4],
                    'Mean TPR': mean_tpr,
                    'FPR': mean_fpr,     
                })
                X = st.session_state.model_results_desc1["X"]
                y = st.session_state.model_results_desc1["y"]
                fig1 = do_plot_roc_curve(tprs_df)
                cols[1].metric("AUROC ± Std", f"{np.mean(aurocs):.3f} ± {std_auroc:.3f}")
                cols[1].success('Model trained!')
                cols[1].altair_chart(fig1, use_container_width=True)

            if 'train_desc2_model_active' not in st.session_state:
                st.session_state['train_desc2_model_active'] = False
            def toggle_train_desc2_model_state():
                st.session_state['train_desc2_model_active'] = not st.session_state['train_desc2_model_active']
            if cols[2].button('🤖 Train a classifier using CC Signatures!', on_click=toggle_train_desc2_model_state):
                if not st.session_state["train_desc2_model_active"]:
                    pass
                else:
                    with st.spinner("Training the model..."):
                        st.session_state.model_results_desc2 = train_acinetobacter_ml_model(st.session_state['desc2_results'].iloc[:, 2:], y)
            if st.session_state["train_desc2_model_active"]:
                if "model_results_desc2" in st.session_state:
                    aurocs = st.session_state.model_results_desc2["aurocs"]
                    std_auroc = np.std(aurocs)
                    tprs = []
                    mean_fpr = np.linspace(0,1,100)
                    for i in st.session_state.model_results_desc2["cv_data"]:
                        fpr, tpr, _ = metrics.roc_curve(i[0],i[1])
                        interp_tpr = np.interp(mean_fpr, fpr, tpr)
                        interp_tpr[0] = 0.0
                        tprs.append(interp_tpr)
                    mean_tpr = np.mean(tprs, axis=0)
                    mean_tpr[-1] = 1.0
                tprs_df = pd.DataFrame({
                    'tpr_cv1': tprs[0],
                    'tpr_cv2': tprs[1],
                    'tpr_cv3': tprs[2],
                    'tpr_cv4': tprs[3],
                    'tpr_cv5': tprs[4],
                    'Mean TPR': mean_tpr,
                    'FPR': mean_fpr,     
                })
                X = st.session_state.model_results_desc2["X"]
                y = st.session_state.model_results_desc2["y"]
                fig1 = do_plot_roc_curve(tprs_df)
                cols[2].metric("AUROC ± Std", f"{np.mean(aurocs):.3f} ± {std_auroc:.3f}")
                cols[2].success('Model trained!')
                cols[2].altair_chart(fig1, use_container_width=True)
                q4f_comb = '  \n'.join(q4_followup)
                cols[0].info(q4f_comb,icon=":material/quiz:")
                st.subheader("Let's dig a bit deeper...")
                cols = st.columns([1,1,1])
                cv_data = st.session_state.model_results_desc2["cv_data"][0]
                cv_data = pd.DataFrame({"ytest":cv_data[0], "ypred": cv_data[1]})
                cols[0].success("Model prediction on test set")
                chart = alt.Chart(cv_data).mark_bar(color="#1D6996").encode(
                    alt.X("ypred", bin=alt.Bin(maxbins=30), axis=alt.Axis(title="Prediction")),
                    y=alt.Y('count()', axis=alt.Axis(title='Counts'))
                    )
                cols[0].altair_chart(chart, use_container_width=True)
                cols[1].success("Threshold of probability")
                proba_cutoff = cols[1].slider("Proba cutoff", min_value=0.1, max_value=1., value=0.5, step=0.001, format="%.2f")
                st.session_state.proba_cutoff = proba_cutoff
                sensitivity, specificity = calculate_sensitivity_recall(cv_data["ytest"], cv_data["ypred"], cutoff = proba_cutoff)
                cols[1].write("Sensitivity: ")
                cols[1].write(sensitivity)
                cols[1].write("Specficity: ")
                cols[1].write(specificity)

                cols[2].success("Contingency table")
                chart = plot_contingency_table(cv_data["ytest"], cv_data["ypred"], cutoff = proba_cutoff)
                cols[2].altair_chart(chart, use_container_width=True)

                if 'final_model' not in st.session_state:
                    st.session_state['final_model'] = False
                def toggle_final_model():
                    st.session_state['final_model'] = not st.session_state['final_model']
                if cols[0].button('💾 Save the model and use it!', on_click=toggle_final_model):
                    pass
                if st.session_state["final_model"]:
                    @st.cache_data
                    def read_library(library_filename):
                        return list(pd.read_csv(os.path.join(data_dir, library_filename))["smiles"])

                    st.divider()
                    st.header("Library selection for prediction")
                    cols = st.columns(5)
                    smiles_list = read_library(library_filenames["Compound library 1"])
                    cols[0].metric("Number of molecules", len(smiles_list))
                    library_molecules_list = [(i, smiles) for i, smiles in enumerate(smiles_list)]
                    num_molecules = len(library_molecules_list)
                    num_chunks_of_4 = (num_molecules + 3) // 4
                    chunk_of_4_index = st.session_state.get('chunk_of_4_index', 0)
                    start_index = chunk_of_4_index * 4
                    end_index = min(start_index + 4, num_molecules)
                    current_chunk_of_4 = library_molecules_list[start_index:end_index]
                    def draw_molecules_in_chunk_of_4(cols, current_chunk):
                        i = 0
                        for m in current_chunk:
                            cols[i+1].image(draw_molecule(m[1]), caption=f"Molecule {m[0]}")
                            i += 1
                    draw_molecules_in_chunk_of_4(cols, current_chunk_of_4)
                    if cols[1].button("View more molecules"):
                        chunk_of_4_index = (chunk_of_4_index + 1) % num_chunks_of_4
                    st.session_state['chunk_of_4_index'] = chunk_of_4_index

                    
                    descriptor_choice = st.radio("Select descriptor for model predictions", ("Morgan", "Chemical Checker"))
                    if 'model_predictions_active' not in st.session_state:
                        st.session_state.dp = None
                        st.session_state['model_predictions_active'] = False
                    def toggle_model_predictions_state():
                        st.session_state['model_predictions_active'] = not st.session_state['model_predictions_active']
                    if st.button(":rocket: Run predictions!", on_click=toggle_model_predictions_state):
                        if not st.session_state["model_predictions_active"]:
                            pass
                        else:
                            with st.spinner("Running model..."):
                                if descriptor_choice == "Morgan":
                                    st.toast("Running the Acinetobacter model")
                                    descs = pd.read_csv(os.path.join("data", "subset250_0_eos4wt0.csv"))
                                    X = descs.iloc[:, 2:]
                                    mdl = st.session_state.model_results_desc1["model"]
                                    abau_preds = predict_acinetobacter_ml_model(X, mdl)
                                elif descriptor_choice == "Chemical Checker":
                                    st.toast("Running the Acinetobacter model")
                                    descs = pd.read_csv(os.path.join("data", "subset250_0_eos4u6p.csv"))
                                    X = descs.iloc[:, 2:]
                                    mdl = st.session_state.model_results_desc2["model"]
                                    abau_preds = predict_acinetobacter_ml_model(X, mdl)
                                abau_df = pd.DataFrame({"smiles": smiles_list, "proba1": abau_preds})
                                cols = st.columns([0.3,0.35,0.35])
                                cols[0].write(abau_df)
                                chart = alt.Chart(abau_df).mark_bar(color="#1D6996").encode(
                                    alt.X("proba1", bin=alt.Bin(maxbins=30), axis=alt.Axis(title="proba1")),
                                    y=alt.Y('count()', axis=alt.Axis(title='Counts'))
                                )
                                cols[1].altair_chart(chart, use_container_width=True)
                                q5_comb = '  \n'.join(q5)
                                cols[2].info(q5_comb, icon=":material/quiz:")