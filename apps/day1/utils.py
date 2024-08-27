import numpy as np
import pandas as pd
import os
from io import BytesIO
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from skimage import img_as_ubyte
import base64
import umap

def clean_df(df):
    """
    df: a Pandas dataframe with a column "SMILES"
    """
    
    #df["CAN_SMILES"] = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) if s is not np.nan else None for s in df["SMILES"]]

    df.drop_duplicates(subset ="SMILES", keep = 'first', inplace = True)
    df = df.dropna()
    df.reset_index(level=None, drop=True, inplace=True, col_level=0, col_fill='')

    return df

def featurize_morgan(smiles_list):
    X = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
            if fp is not None:
                X.append(fp)
            else:
                X.append(None)
        else:
            X.append(None)
    return X

def create_umap(fp_list):
    umap_transformer = umap.UMAP(a=0.001, b=1.5, min_dist=0.1)
    return umap_transformer.fit_transform(fp_list)

def image_formatter(smiles):
    mol = Chem.MolFromSmiles(smiles)
    image = Draw.MolToImage(mol, size=(150,150))
    with BytesIO() as buffer:
        image.save(buffer, 'png')
        data = base64.encodebytes(buffer.getvalue()).decode('utf-8')
    return f"data:image/png;base64,{data}"
