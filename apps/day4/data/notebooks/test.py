from rdkit import Chem
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
from rdkit.Chem import AllChem
import numpy as np

df = pd.read_csv("provisional.csv")

pos_smiles = df[df["pchembl"] >= 5]["smiles"].tolist()
pos_smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in pos_smiles]

neg_smiles = pd.read_csv("../raw/drugbank_smiles.csv")["Smiles"].tolist()
neg_smiles = [x for x in neg_smiles if x not in pos_smiles]
neg_smiles_ = []
for smi in neg_smiles:
    try:
        x = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    except:
        continue
    neg_smiles_ += [x]
print(len(neg_smiles_))
neg_smiles = neg_smiles_

smiles_list = pos_smiles + neg_smiles
y = [1] * len(pos_smiles) + [0] * len(neg_smiles)
idxs = [i for i in range(len(y))]
random.shuffle(idxs)
smiles_list = [smiles_list[i] for i in idxs]
y = [y[i] for i in idxs]


def fingerprints(smiles_list):

    molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    # Define Morgan fingerprint parameters
    radius = 2
    nBits = 1024

    # Calculate fingerprints
    fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits) for mol in molecules]

    # Convert fingerprints to numpy array
    fingerprint_matrix = np.array([np.array(fingerprint) for fingerprint in fingerprints])

    return fingerprint_matrix

X = fingerprints(smiles_list)
model = RandomForestClassifier()
model.fit(X, y)

library_smiles = pd.read_csv("../processed/all_mols.csv")["smiles"].tolist()
X = fingerprints(library_smiles)
y_hat = model.predict_proba(X)[:,1]

print(y_hat)
print(np.max(y_hat))