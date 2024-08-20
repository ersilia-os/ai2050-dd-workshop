import requests
import os
import csv
import random

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