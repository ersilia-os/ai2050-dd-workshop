import sys
import os
import subprocess
import tempfile
import pandas as pd
import json
import csv
import shutil
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

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
    return os.path.join(root, "..", "cache", f"sampling_{input_inchikey}_{model_id}_{round}.json")

def get_prediction_csv_filename(model_id):
    return os.path.join(root, "..", "cache", f"prediction_{input_inchikey}_{model_id}.csv")

def ersilia_sampling_runner(model_id, input_file, round=0):
    print("Running model:", model_id, "Round:", round)
    output_file = os.path.join(tmp_dir, f"output_{model_id}_{round}.csv")
    individual_script_file = os.path.join(tmp_dir, "individual_script.sh")
    with open(individual_script_file, "w") as f:
        f.write("ersilia -v serve {0}\n".format(model_id))
        f.write("ersilia -v run -i {0} -o {1}\n".format(input_file, output_file))
        f.write("ersilia close\n")
    cmd = f"bash {individual_script_file}"
    print("Running command:", cmd)
    subprocess.run(cmd, shell=True)
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
        f.write("ersilia -v serve {0}\n".format(model_id))
        f.write("ersilia -v run -i {0} -o {1}\n".format(internal_input_file, output_file))
        f.write("ersilia close\n")
    cmd = f"bash {individual_script_file}"
    print("Running command:", cmd)
    subprocess.run(cmd, shell=True)
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
smiles_list = [input_smiles]
for model_id in sampling_models:
    for i in range(3):
        json_file = get_sampling_json_filename(model_id, i)
        if not os.path.exists(json_file):
            continue
        with open(json_file, "r") as f:
            smiles_list += json.load(f)
        print(len(smiles_list))

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

print("All done!")
print("Cleaning up")
shutil.rmtree(tmp_dir)

print("Assembling output in dataframes")
cache_folder = os.path.abspath(os.path.join(root, "..", "cache"))

print("Starting with the sampling results")
R = []
for fn in os.listdir(cache_folder):
    if not fn.startswith("sampling"):
        continue
    with open(os.path.join(cache_folder, fn), "r") as f:
        data = json.load(f)
        _, input_inchikey, model_id, round = fn.split(".json")[0].split("_")
        round = int(round)
        for smiles in data:
            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
            R += [[model_id, round, inchikey, smiles]]

df = pd.DataFrame(R, columns = ["model_id", "round", "inchikey", "smiles"])
file_name = os.path.join(root, "..", "cache_lite", f"sampling_{input_inchikey}.csv")
if os.path.exists(file_name):
    df_0 = pd.read_csv(file_name)
    df = pd.concat([df_0, df]).reset_index(drop=True)
df = df.drop_duplicates(inplace=False)
df.to_csv(file_name, index=False)

print("Finishing with the predition results")
dfs = []
model_molecule = []
for fn in os.listdir(cache_folder):
    if not fn.startswith("prediction"):
        continue
    df = pd.read_csv(os.path.join(cache_folder, fn))
    _, input_inchikey, model_id = fn.split(".csv")[0].split("_")
    dfs += [df]
    model_molecule += [(model_id, input_inchikey)]

inchikeys = dfs[0]["key"].tolist()
smiles = dfs[0]["input"].tolist()
for df in dfs:
    assert inchikeys == df["key"].tolist(), "Fatal error! Molecules are not the same between models."

df = pd.DataFrame({"inchikey": inchikeys, "smiles": smiles})

for i, df_ in enumerate(dfs):
    print(df)
    columns = list(df_.columns)[2:]
    df_ = df_[columns]
    model_id = model_molecule[i][0]
    rename = dict((c, model_id + "_" + c) for c in columns)
    df_ = df_.rename(columns=rename, inplace=False)
    df = pd.concat([df, df_], axis=1)
    columns = ["inchikey", "smiles"]
    columns += [x for x in columns if x in activity_models]
    columns += [x for x in columns if x in adme_models]
    df = df[columns]

df.to_csv(os.path.join(root, "..", "cache_lite", f"predictions_{input_inchikey}.csv"), index=False)
    