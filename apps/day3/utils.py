import requests
import os
import csv
import random
import pandas as pd
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

root = os.path.dirname(os.path.abspath(__file__))


def get_model_info_from_github(model_id):
    url = "https://raw.githubusercontent.com/ersilia-os/{0}/main/metadata.json".format(model_id)
    response = requests.get(url)
    return response.json()


def random_molecules(n=100):
    smiles_list = []
    with open(os.path.join(root, "data", "example.csv")) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            smiles_list.append(row[0])
    return random.choices(smiles_list, k=n)


def basic_molecules_dataframe(smiles_list, ref_smiles):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    similarity_scores = []
    fpgen = AllChem.GetRDKitFPGenerator()
    ref_fp = fpgen.GetFingerprint(ref_mol)
    for mol in mols:
        mol_fp = fpgen.GetFingerprint(mol)
        similarity = DataStructs.FingerprintSimilarity(ref_fp, mol_fp)
        similarity_scores.append(round(similarity, 3))
    inchikeys = [Chem.MolToInchiKey(mol) for mol in mols]
    df = pd.DataFrame({"InChIKey": inchikeys})
    df["SMILES"] = smiles_list
    df["Tanimoto Coeff"] = similarity_scores
    df["MolWeight"] = [round(Descriptors.MolWt(mol), 1) for mol in mols]
    df["LogP"] = [round(Descriptors.MolLogP(mol), 3) for mol in mols]
    df["QED"] = [round(Chem.QED.qed(mol), 3) for mol in mols]
    return df.sort_values("Tanimoto Coeff", ascending=False).reset_index(drop=True)