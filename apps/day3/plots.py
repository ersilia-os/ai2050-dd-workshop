from rdkit import Chem
from rdkit.Chem import Draw
import mols2grid
import plotly.express as px
import streamlit as st


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
    markdown_text += f"- [**Publication**]({model_info['Publication']})\n"
    markdown_text += f"- [**Source Code**]({model_info['Source Code']})\n"
    github_repo = f"https://github.com/ersilia-os/{model_id}"
    markdown_text += f"- [**Ersilia Repository**]({github_repo})" + "\n"
    markdown_text += "\n"
    markdown_text += model_info['Description'] + "\n"
    container.markdown(markdown_text)

def plot_single_molecule(container, smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    container.image(img)

def basic_molecule_card(container, smiles, df):
    plot_single_molecule(container, smiles)
    df = df[df["SMILES"] == smiles]
    if df.shape[0] == 0:
        return
    text = "- **InChIKey**: {0}\n".format(df["InChIKey"].values[0])
    text += "- **SMILES**: {0}\n".format(df["SMILES"].values[0])
    text += "- **Tanimoto Coeff**: {0}\n".format(df["Tanimoto Coeff"].values[0])
    text += "- **MolWeight**: {0}\n".format(df["MolWeight"].values[0])
    text += "- **LogP**: {0}\n".format(df["LogP"].values[0])
    text += "- **QED**: {0}\n".format(df["QED"].values[0])
    container.markdown(text)

def molecule_card(container, smiles, df, activity_column=None):
    plot_single_molecule(container, smiles)
    df = df[df["SMILES"] == smiles]
    if df.shape[0] == 0:
        return
    text = "- **InChIKey**: {0}\n".format(df["InChIKey"].values[0])
    text += "- **SMILES**: {0}\n".format(df["SMILES"].values[0])
    text += "- **Tanimoto Coeff**: {0}\n".format(df["Tanimoto Coeff"].values[0])
    if activity_column is None:
        text += "- **Activity**: {0}\n".format(df["Activity"].values[0])
    else:
        text += "- **{0}**: {1:.3f}\n".format(activity_column, df[activity_column].values[0])
    cols = list(df.columns)[4:]
    for col in cols:
        text += "- **{0}**: {1}\n".format(col, df[col].values[0])
    container.markdown(text)

def basic_mols2grid_plot(df):
    html_str = mols2grid.display(df, selection=False, template="interactive")._repr_html_()
    mols2grid.get_selection()
    st.components.v1.html(html_str, height=350, scrolling=True)

def histogram(container, df, column):
    fig = px.histogram(df, x=column)
    container.plotly_chart(fig, theme="streamlit")
