import requests
import os
import csv
import random
import collections
import pandas as pd
import numpy as np
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

root = os.path.dirname(os.path.abspath(__file__))
cache_dir = os.path.abspath(os.path.join(root, "data", "cache_lite"))

deterministic_models = [
    "eos9ueu",
    "eos1d7r",
    "eos3kcw",
]

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


def cached_molecules(inchikey, model_id, n=100):
    cache_file = os.path.join(cache_dir, "sampling_{0}.csv".format(inchikey))
    if model_id in deterministic_models:
        if not os.path.exists(cache_file):
            return None
        print(cache_file)
        df = pd.read_csv(cache_file)
        df = df[df["model_id"] == model_id].reset_index(drop=True)
        print(df)
        ik2smi = {}
        for v in df[["inchikey", "smiles"]].values:
            ik2smi[v[0].split("-")[0]] = v[1]
        rounds = collections.defaultdict(list)
        round_ids = sorted(set(df["round"]))
        for round in round_ids:
            rounds[round] = [x.split("-")[0] for x in df[df["round"] == round]["inchikey"].tolist()]
        sampled_molecules = []
        sampled_molecules_set = set()
        for _ in range(1000000):
            for k, v in rounds.items():
                if len(v) > 0:
                    x = v.pop()
                    if x not in sampled_molecules_set:
                        sampled_molecules.append(ik2smi[x])
                        sampled_molecules_set.add(x)
                        if len(sampled_molecules) >= n:
                            return sampled_molecules
        return sampled_molecules
    else:
        if not os.path.exists(cache_file):
            return None
        df = pd.read_csv(cache_file)
        df = df[df["model_id"] == model_id].reset_index(drop=True)
        smiles_list = set(df["smiles"].tolist())
        return random.sample(smiles_list, k=np.min([n, len(smiles_list)]))


def random_activity(n):
    return [round(random.uniform(0, 1), 2) for _ in range(n)]


def random_adme_percentiles(n, m):
    R = []
    for _ in range(n):
        r = []
        for _ in range(m):
            r.append(round(random.uniform(0, 1)*100, 2))
        R.append(r)
    return R


def cached_predictions(inchikey, model_id, smiles_list):
    cache_file = os.path.join(cache_dir, "predictions_{0}.csv".format(inchikey))
    if not os.path.exists(cache_file):
        return None
    df = pd.read_csv(cache_file)
    all_columns = list(df.columns)
    columns = ["inchikey", "smiles"] + [x for x in all_columns if x.startswith(model_id)]
    df = df[columns]
    data = {}
    for r in df.values:
        data[r[1]] = r
    R = []
    for smiles in smiles_list:
        R += [data[smiles]]
    df = pd.DataFrame(R, columns=columns)
    return df


def smiles_to_inchikey(smiles):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))


def filter_smiles_library(inchikey, smiles_list):
    cache_file = os.path.join(cache_dir, "predictions_{0}.csv".format(inchikey))
    if not os.path.exists(cache_file):
        return None
    df = pd.read_csv(cache_file)
    all_smiles = set(df["smiles"].tolist())
    smiles_list = [x for x in smiles_list if x in all_smiles]
    return smiles_list


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