from rdkit import Chem
from rdkit.Chem import Draw

def markdown_card(container, model_id, models_info):
    model_info = models_info[model_id]
    markdown_text = f"### :robot_face: {model_info['Title']}\n"
    markdown_text += f"- **EOS ID**: {model_id}\n"
    markdown_text += f"- **Slug**: {model_info['Slug']}\n"
    markdown_text += f"- **Mode**: {model_info['Mode']}\n"
    task = ", ".join(model_info['Task'])
    markdown_text += f"- **Task**: {task}\n"
    inp = ", ".join(model_info['Input'])
    markdown_text += f"- **Input**: {inp}\n"
    out = ", ".join(model_info['Output'])
    markdown_text += f"- **Output**: {out}\n"
    markdown_text += f"- **Publication**: {model_info['Publication']}\n"
    markdown_text += f"- **Source Code**: {model_info['Source Code']}\n"
    github_repo = f"https://github.com/ersilia-os/{model_id}"
    markdown_text += f"- **Ersilia Repository**: {github_repo}" + "\n"
    markdown_text += "\n"
    markdown_text += model_info['Description'] + "\n"
    container.markdown(markdown_text)


def plot_single_molecule(container, smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    container.image(img)