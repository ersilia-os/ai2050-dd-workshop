import os
import numpy as np
import pandas as pd
from lol import LOL


def load_acinetobacter_training_data(datapath):
    df = pd.read_csv(datapath)
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    df.reset_index(drop=True, inplace=True)
    df["Mean"] = pd.to_numeric(df["Mean"])
    return df

def binarize_acinetobacter_data(data, cutoff):
    y = [1 if i <= cutoff else 0 for i in data["Mean"]]
    columns = list(data.columns)[:5]
    data = data[columns]
    data["Binary"] = y
    return data

def lolp_reducer(X, y):
    reducer = LOL(100)
    X_ = reducer.fit_transform(X, y)
    results = {
        "reducer": reducer,
        "X": X_,
        "y": y
    }
    return results