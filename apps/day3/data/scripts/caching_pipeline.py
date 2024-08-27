import sys
import os
import subprocess
import tempfile
import pandas as pd
from rdkit import Chem
import json
import csv

input_smiles = sys.argv[1]
input_inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(input_smiles))

tmp_dir = tempfile.mkdtemp()
print("Working in directory:", tmp_dir)
print("Input SMILES:", input_smiles)

root = os.path.dirname(os.path.abspath(__file__))

sampling_models = [
    "eos8fma",
    "eos4q1a",
    "eos6ost",
    "eos9ueu",
    "eos1d7r",
    "eos3kcw",
    "eos694w",
]

activity_models = [
    "eos5xng",
    "eos4e40",
    "eos9f6t",
]

adme_models = [
    "eos7d58"
]

print("Creating input file")
input_file = os.path.join(tmp_dir, "input.csv")

with open(input_file, "w") as f:
    f.write("smiles\n")
    f.write(input_smiles)

def get_sampling_json_filename(model_id, round):
    return os.path.join(root, "..", "data", f"sampling_{input_inchikey}_{model_id}_{round}.json")

def get_prediction_csv_filename(model_id):
    return os.path.join(root, "..", "data", f"prediction_{input_inchikey}_{model_id}.csv")

def ersilia_sampling_runner(model_id, input_file, round=0):
    print("Running model:", model_id, "Round:", round)
    output_file = os.path.join(tmp_dir, f"output_{model_id}_{round}.csv")
    individual_script_file = os.path.join(tmp_dir, "individual_script.sh")
    with open(individual_script_file, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("conda activate ersilia\n")
        f.write("ersilia serve {0}\n".format(model_id))
        f.write("ersilia run -i {0} -o {1}\n".format(input_file, output_file))
        f.write("ersilia close\n")
    cmd = f"bash {individual_script_file}"
    print("Running command:", cmd)
    subprocess.run(cmd, shell=True).wait()
    print("Output file:", output_file)
    return output_file

def ersilia_prediction_runner(model_id, smiles_list):
    internal_input_file = os.path.join(tmp_dir, "internal_input.csv")
    with open(internal_input_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles"])
        for smiles in smiles_list:
            writer.writerow([smiles])
    print("Running model:", model_id)
    output_file = os.path.join(tmp_dir, f"output_{model_id}.csv")
    individual_script_file = os.path.join(tmp_dir, "individual_script.sh")
    with open(individual_script_file, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("conda activate ersilia\n")
        f.write("ersilia serve {0}\n".format(model_id))
        f.write("ersilia run -i {0} -o {1}\n".format(input_file, output_file))
        f.write("ersilia close\n")
    cmd = f"bash {individual_script_file}"
    print("Running command:", cmd)
    subprocess.run(cmd, shell=True).wait()
    print("Output file:", output_file)
    return output_file

def convert_tmp_sampling_output_to_json(output_file, model_id, round):
    json_file = get_sampling_json_filename(model_id, round)
    print("Converting output to JSON:", output_file)
    df = pd.read_csv(output_file)
    for v in df.values:
        v = v[2:]
    smiles_list = []
    for x in v:
        x = str(x)
        mol = Chem.MolFromSmiles(x)
        if mol is None:
            continue
        smiles_list.append(x)
    with open(json_file, "w") as f:
        json.dump(smiles_list, f, indent=4)

def keep_tmp_prediction_output(output_file, model_id):
    csv_file = get_prediction_csv_filename(model_id)
    print("Keeping output in:", output_file)
    df = pd.read_csv(output_file)
    df.to_csv(csv_file, index=False)

print("Running sampling models")
for model_id in sampling_models:
    print("Model:", model_id)
    for i in range(3):
        print("Round:", i)
        json_file = get_sampling_json_filename(model_id, i)
        if os.path.exists(json_file):
            continue
        output_file = ersilia_sampling_runner(model_id, input_file, i)
        convert_tmp_sampling_output_to_json(output_file, model_id, i)

print("Collecting sampling results")
for model_id in sampling_models:
    for i in range(3):
        json_file = get_sampling_json_filename(model_id, i)
        if not os.path.exists(json_file):
            continue
        with open(json_file, "r") as f:
            smiles_list = json.load(f)

print("Running activity models")
for model_id in activity_models:
    print("Model:", model_id)
    csv_file = get_prediction_csv_filename(model_id)
    if os.path.exists(csv_file):
        print("Already exists. Skipping.")
        continue
    output_file = ersilia_prediction_runner(model_id, smiles_list)
    keep_tmp_prediction_output(output_file, model_id)

print("Running ADME models")
for model_id in adme_models:
    print("Model:", model_id)
    csv_file = get_prediction_csv_filename(model_id)
    if os.path.exists(csv_file):
        print("Already exists. Skipping.")
        continue
    output_file = ersilia_prediction_runner(model_id, smiles_list)
    keep_tmp_prediction_output(output_file, model_id)