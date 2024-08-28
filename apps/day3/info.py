import os
import collections
import pandas as pd

root = os.path.abspath(os.path.dirname(__file__))

about = [
    "This app is part of the Schmidt Futures AI2050 program.",
    "The workshop has been jointly developed by the [Ersilia Open Source Initiative](https://ersilia.io) and the [H3D Foundation](https://h3dfoundation.org/).",
    "For more information about this workshop, please see code and data in this [GitHub repository](https://github.com/ersilia-os/ai2050-dd-workshop).",
    "The app is intended for educational purposes only. If you have a more advanced dataset in mind or a use case for your research, please contact us at: [hello@ersilia.io](mailto:hello@ersilia.io)."
]

intro = """
In this session, we will explore sampling methods to expand the chemical space around seed (hit) molecules in a _Burkholderia cenocepacia_ screening assay. We will try and critically assess a few exemplary methods, both for generative chemistry and for querying (ultra) large scale chemical libraries.
""".strip()

step_1_explanation = """
1. **Seed molecules**: Select a seed molecule from the list of seed molecules.
1. **Ersilia Model Hub identifiers**: Select an Ersilia model from the list of models. Read model descriptions to understand the model's capabilities.
1. **Sample molecules!**: Click this button to sample molecules from the selected model.
1. **Candidates library**: Create a small chemical library by **copy-pasting** SMILES strings from the table on the left into the input text area below.
""".strip()

step_1_questions = """
- What models give deterministic results? What models give stochastic results?
- Qualitatively, how diverse are the generated molecules?
- How do the sampled molecules compare to the seed molecule?
- Do model results overlap?
- What is the rationale for selecting molecules for the candidates library?
""".strip()

step_2_explanation = """
1. **Select a primary activity predictor**: Select the activity prediction model that is best suited to your task. Read the model description carefully.
1. **Select ADMET properties**: Select up to 5 ADMET properties to predict.
1. **Predict properties**: Click the button to predict the properties of the molecules in your wishlist. Results will be visible in a table below.
""".strip()

step_2_questions = """
- What is the rationale for selecting the activity predictor?
- What is the rationale for selecting the ADMET properties?
- Are there other properties that you would like to predict?
""".strip()

step_3_explanation = """
1. **Explore molecules**: Explore individual molecules by copy-pasting their SMILES strings into the input text area.
1. **Inspect property distributions**: Critically analyze the distribution of properties one by one.
1. **Make the final selection**: Select 5 candidates for further analysis. Discuss this selection with your colleagues.
""".strip()

step_3_questions = """
- What are the most important properties for your task?
- How do ADMET properties of the molecules in your wishlist compare to the properties of known drugs?
- How did you select your top 5 candidates?
- Was it easy to reach a consensus with your colleagues?
- What are the limitations of the models used in this session?
- What would you do next?
""".strip()

adme_warning_message = """
ADME properties will be given as percentiles with respect to approved drugs in DrugBank. A value of 100 means that the property is higher than all approved drugs in DrugBank. A value of 0 means that the property is lower than all approved drugs in DrugBank. A value of 50 means that the property is in the middle of the distribution of approved drugs in DrugBank.
""".strip()

model_urls = {
    "eos8fma": None,
    "eos4q1a": None,
    "eos6ost": None,
    "eos9ueu": None,
    "eos1d7r": None,
    "eos3kcw": None,
    "eos694w": None,
}

activity_models_urls = {
    "eos5xng": None,
    "eos4e40": None,
    "eos9f6t": None,
}

adme_models_urls = {
    "eos7d58": None
}

deterministic_models = [
    "eos9ueu",
    "eos1d7r",
    "eos3kcw",
]

ds = pd.read_csv(os.path.join(root, "data", "atb_screening_data_top100.csv")).head(4)
seed_molecules = collections.OrderedDict()
for v in ds.values:
    seed_molecules[v[0]] = v[1]

raw_adme_header = "molecular_weight,logP,hydrogen_bond_acceptors,hydrogen_bond_donors,Lipinski,QED,stereo_centers,tpsa,AMES,BBB_Martins,Bioavailability_Ma,CYP1A2_Veith,CYP2C19_Veith,CYP2C9_Substrate_CarbonMangels,CYP2C9_Veith,CYP2D6_Substrate_CarbonMangels,CYP2D6_Veith,CYP3A4_Substrate_CarbonMangels,CYP3A4_Veith,Carcinogens_Lagunin,ClinTox,DILI,HIA_Hou,NR-AR-LBD,NR-AR,NR-AhR,NR-Aromatase,NR-ER-LBD,NR-ER,NR-PPAR-gamma,PAMPA_NCATS,Pgp_Broccatelli,SR-ARE,SR-ATAD5,SR-HSE,SR-MMP,SR-p53,Skin_Reaction,hERG,Caco2_Wang,Clearance_Hepatocyte_AZ,Clearance_Microsome_AZ,Half_Life_Obach,HydrationFreeEnergy_FreeSolv,LD50_Zhu,Lipophilicity_AstraZeneca,PPBR_AZ,Solubility_AqSolDB,VDss_Lombardo,molecular_weight_drugbank_approved_percentile,logP_drugbank_approved_percentile,hydrogen_bond_acceptors_drugbank_approved_percentile,hydrogen_bond_donors_drugbank_approved_percentile,Lipinski_drugbank_approved_percentile,QED_drugbank_approved_percentile,stereo_centers_drugbank_approved_percentile,tpsa_drugbank_approved_percentile,AMES_drugbank_approved_percentile,BBB_Martins_drugbank_approved_percentile,Bioavailability_Ma_drugbank_approved_percentile,CYP1A2_Veith_drugbank_approved_percentile,CYP2C19_Veith_drugbank_approved_percentile,CYP2C9_Substrate_CarbonMangels_drugbank_approved_percentile,CYP2C9_Veith_drugbank_approved_percentile,CYP2D6_Substrate_CarbonMangels_drugbank_approved_percentile,CYP2D6_Veith_drugbank_approved_percentile,CYP3A4_Substrate_CarbonMangels_drugbank_approved_percentile,CYP3A4_Veith_drugbank_approved_percentile,Carcinogens_Lagunin_drugbank_approved_percentile,ClinTox_drugbank_approved_percentile,DILI_drugbank_approved_percentile,HIA_Hou_drugbank_approved_percentile,NR-AR-LBD_drugbank_approved_percentile,NR-AR_drugbank_approved_percentile,NR-AhR_drugbank_approved_percentile,NR-Aromatase_drugbank_approved_percentile,NR-ER-LBD_drugbank_approved_percentile,NR-ER_drugbank_approved_percentile,NR-PPAR-gamma_drugbank_approved_percentile,PAMPA_NCATS_drugbank_approved_percentile,Pgp_Broccatelli_drugbank_approved_percentile,SR-ARE_drugbank_approved_percentile,SR-ATAD5_drugbank_approved_percentile,SR-HSE_drugbank_approved_percentile,SR-MMP_drugbank_approved_percentile,SR-p53_drugbank_approved_percentile,Skin_Reaction_drugbank_approved_percentile,hERG_drugbank_approved_percentile,Caco2_Wang_drugbank_approved_percentile,Clearance_Hepatocyte_AZ_drugbank_approved_percentile,Clearance_Microsome_AZ_drugbank_approved_percentile,Half_Life_Obach_drugbank_approved_percentile,HydrationFreeEnergy_FreeSolv_drugbank_approved_percentile,LD50_Zhu_drugbank_approved_percentile,Lipophilicity_AstraZeneca_drugbank_approved_percentile,PPBR_AZ_drugbank_approved_percentile,Solubility_AqSolDB_drugbank_approved_percentile,VDss_Lombardo_drugbank_approved_percentile"
raw_adme_header = raw_adme_header.split(",")

db_adme_header = []
for h in raw_adme_header:
    if "drugbank_approved_percentile" not in h:
        continue
    else:
        db_adme_header.append(h)

adme_col_new2old = collections.OrderedDict()
for old in db_adme_header:
    new = old.replace("_drugbank_approved_percentile", "")
    new = new.replace("_", " ")
    new = new[0].upper() + new[1:]
    adme_col_new2old[new] = old
