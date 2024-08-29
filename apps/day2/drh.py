import pandas as pd
from rdkit import Chem

# Step 1: Load the data, skipping header lines that start with '!'
file_path = "data/repurposing_samples_20240610.txt"
df = pd.read_csv(file_path, sep='\t', comment='!')

# Step 2: Select the columns 'smiles', 'InChIKey', and 'pubchem_cid'
df_filtered = df[['smiles', 'InChIKey', 'pubchem_cid']]

# Step 3: Save the filtered data to a CSV file
output_path = "data/drugrepurposinghub.csv"
df_filtered = df_filtered.drop_duplicates(keep = "first")
def is_valid_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None
df_filtered = df_filtered[df_filtered['smiles'].apply(is_valid_smiles)]
df_filtered.to_csv(output_path, index=False)

print(f"Data saved to {output_path}")
