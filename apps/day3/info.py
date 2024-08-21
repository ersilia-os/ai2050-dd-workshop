import collections

about = [
    "This app is part of the Schmidt Futures AI2050 program.",
    "The workshop been jointly developed by the [Ersilia Open Source Initiative](https://ersilia.io) and the [H3D Foundation](https://h3dfoundation.org/).",
    "For more information about this workshop, please see code and data in this [GitHub repository](https://github.com/ersilia-os/ai2050-dd-workshop).",
    "If you have a more advanced dataset in mind or a use case for your research, please contact us at: [hello@ersilia.io](mailto:hello@ersilia.io)."
]

intro = """
In this workshop, we will explore sampling methods to expand the chemical space around a seed (hit) molecule. We will try and critically assess a few exemplary methods, both for generative chemistry and for querying (ultra) large scale chemical libraries.
""".strip()

step_1_explanation = """
- **Seed molecules**: Select a seed molecule from the list of seed molecules.
- **Ersilia Model Hub identifiers**: Select a model from the list of models.
- **Sample molecules!**: Click this button to sample molecules from the selected model.
""".strip()

step_2_explanation = """
- **Seed molecules**: Select a seed molecule from the list of seed molecules.
- **Ersilia Model Hub identifiers**: Select a model from the list of models.
- **Sample molecules!**: Click this button to sample molecules from the selected model.
""".strip()

step_3_explanation = """
- **Seed molecules**: Select a seed molecule from the list of seed molecules.
- **Ersilia Model Hub identifiers**: Select a model from the list of models.
- **Sample molecules!**: Click this button to sample molecules from the selected model.
""".strip()

adme_warning_message = """
ADME properties will be given as percentiles with respect to approved drugs in DrugBank. A value of 100 means that the property is higher than all approved drugs in DrugBank. A value of 0 means that the property is lower than all approved drugs in DrugBank. A value of 50 means that the property is in the middle of the distribution of approved drugs in DrugBank."
""".strip()

model_urls = {
    "eos8fma": None,
    "eos4q1a": None,
    "eos6ost": None,
    "eos9ueu": None,
    "eos1d7r": None,
    "eos3kcw": None,
}

activity_models_urls = {
    "eos4e40": None,
    "eos4e41": None,
}

adme_models_urls = {
    "eos7d58": None
}

seed_molecules = {
    "Molecule 1": "CC1=CC(=O)C(=C(C1=O)O)C",
    "Molecule 2": "C1=CC=C(C=C1)C(=O)O",
    "Molecule 3": "CC(C)C1=CC(=C(C=C1)O)C",
    "Molecule 4": "C1=CC=C(C=C1)C(=O)O",
}

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
