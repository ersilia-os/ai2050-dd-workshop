import pandas as pd
from rdkit import Chem
from standardiser import standardise

file_path = "data/repurposing_samples_20240610.txt"
df = pd.read_csv(file_path, sep='\t', comment='!')
df_ = df[['smiles', 'InChIKey', 'pubchem_cid']]
print(df_.shape)

output_path = "data/drugrepurposinghub.csv"
df_ = df_.drop_duplicates(keep = "first")
print(df_.shape)

smiles = []
for s in df_["smiles"]:
    mol = Chem.MolFromSmiles(s)
    if mol is not None:
        try:
            mol_ = standardise.run(mol)
            s_ = Chem.MolToSmiles(mol_)
        except:
            print("Non standard:", s)
            s_ = None
    else:
        print("No rdkit:" ,s)
        s_ = None
    smiles += [s_]

df_["smiles_ok"] = smiles
df_ok = df_[~df_["smiles_ok"].isna()]
print(df_ok.shape)

df_ok.to_csv(output_path, index=False)

