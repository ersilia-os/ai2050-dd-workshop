import streamlit as st
import pandas as pd
import csv
import io
import requests
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import matplotlib.pyplot as plt
from xlsxwriter.utility import xl_col_to_name

st.set_page_config(
    page_title="Assemble model predictions from Ersilia",
    page_icon=":molecule:",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.title("Assemble model predictions from Ersilia")

uploaded_input = st.file_uploader(
    "Upload a list of SMILES strings",
    accept_multiple_files=False,
    help="You can upload a list of SMILES strings",
)

def clean_dataframe(model_id, df):
    columns = list(df.columns)
    columns_to_remove = []
    for col in columns:
        if col.lower() == "input":
            columns_to_remove.append(col)
        if col.lower() == "smiles":
            columns_to_remove.append(col)
        if col.lower() == "key":
            columns_to_remove.append(col)
        if col.lower() == "inchikey":
            columns_to_remove.append(col)
        if "_drugbank_" in col.lower():
            columns_to_remove.append(col)
    if columns_to_remove:
        df = df.drop(columns=columns_to_remove)
    columns = list(df.columns)
    columns = [col + "_" + model_id for col in columns]
    df = df.rename(columns=dict(zip(list(df.columns), columns)))
    return df

def get_slug(model_id):
    file_name = "https://raw.githubusercontent.com/ersilia-os/{0}/refs/heads/main/README.md".format(model_id)
    response = requests.get(file_name)
    if response.status_code == 200:
        file_text = response.text
    else:
        file_text = ""
    for line in file_text.split("\n"):
        if line.startswith("- **Slug:**"):
            slug = line.split("`")[1]
            return slug
    st.error(f"Could not find slug for model {model_id}. Please check the model ID.")
        
def get_columns(model_id):
    csv_file = "https://raw.githubusercontent.com/ersilia-os/{0}/refs/heads/main/model/framework/columns/run_columns.csv".format(model_id)
    df = pd.read_csv(csv_file)
    return df

smiles_list = []
if uploaded_input is not None:
    smiles_list = []
    text_io = io.StringIO(uploaded_input.getvalue().decode("utf-8"))
    reader = csv.reader(text_io)

    for row in reader:
        if row:
            smiles_list += [row[0]]

if len(smiles_list) > 0:
    
    uploaded_files = st.file_uploader(
        "Upload up to 10 files produced by Ersilia", 
        type=["csv"], 
        accept_multiple_files=True, 
        help="You can upload up to 5 files."
    )

    if uploaded_files:
        if len(uploaded_files) > 10:
            st.error("Please upload no more than 10 files.")
        else:
            st.success(f"{len(uploaded_files)} file(s) uploaded successfully.")
            model_ids = []
            dataframes = []
            for file in uploaded_files:
                fn = file.name
                if "eos" not in fn or not fn.endswith(".csv"):
                    st.error(f"Invalid file name: {fn}. Please upload files with names containing 'eos' identifier and ending with '.csv'.")
                    st.stop()
                for x in fn.split("_"):
                    if x.startswith("eos") and len(x) == 7:
                        model_id = x
                model_ids += [model_id]
                df = pd.read_csv(file)
                dataframes += [df]

            if len(model_ids) != len(set(model_ids)):
                st.error("Duplicate model IDs found in the uploaded files. Please ensure each file has a unique model ID.")
                st.stop()

            data = {
                "smiles": smiles_list
            }
            data = pd.DataFrame(data)

            for model_id, df in zip(model_ids, dataframes):
                df = clean_dataframe(model_id, df)
                data = pd.concat([data, df], axis=1)

            columns = list(data.columns)[1:]

            df = data.copy()

            cpd_ids = ["mol-{0}".format(str(i).zfill(4)) for i in range(len(df))]

            df["identifier"] = cpd_ids

            # Convert SMILES to RDKit Mol objects
            df["molecule"] = df["smiles"].apply(lambda x: Chem.MolFromSmiles(x))

            df = df[["identifier", "smiles", "molecule"] + columns]

            # Create Excel file with embedded images
            output = BytesIO()
            writer = pd.ExcelWriter(output, engine="xlsxwriter")
            df.to_excel(writer, index=False, startrow=0, sheet_name="Ersilia Results")
            workbook = writer.book
            worksheet = writer.sheets["Ersilia Results"]

            # Set first column (A) to be wide and left aligned
            smiles_col_width = 40
            left_align_format = workbook.add_format({'align': 'left', 'valign': 'vcenter', 'text_wrap': True})

            # Apply left_align_format to all columns except C (molecule image)
            
            num_cols = len(df.columns)
            for col_idx in range(num_cols):
                col_letter = xl_col_to_name(col_idx)
                if col_letter != 'C':
                    worksheet.set_column(f'{col_letter}:{col_letter}', None, left_align_format)

            worksheet.set_column('B:B', smiles_col_width, left_align_format)

            # Set column width and row height to fit images
            img_width, img_height = 100, 100  # pixels
            worksheet.set_column('C:C', img_width / 7.15)  # approx conversion

            for idx, mol in enumerate(df["molecule"], start=1):
                worksheet.set_row(idx, img_height*0.8)  # approx conversion
                if mol:
                    img = Draw.MolToImage(mol, size=(img_width, img_height))
                    img_io = BytesIO()
                    img.save(img_io, format='PNG')
                    img_io.seek(0)
                    # Center image in cell
                    x_offset = 3
                    y_offset = 3
                    worksheet.insert_image(
                        f'C{idx+1}',
                        'molecule.png',
                        {
                            'image_data': img_io,
                            'x_offset': max(x_offset, 0),
                            'y_offset': max(y_offset, 0),
                            'x_scale': 1.0,
                            'y_scale': 1.0,
                        }
                    )

            # Generate a color palette for the model_ids
            palette = plt.get_cmap('tab10')
            model_colors = {model_id: palette(i % 10) for i, model_id in enumerate(model_ids)}
            # Convert RGBA to hex
            for k, v in model_colors.items():
                model_colors[k] = '#%02x%02x%02x' % tuple(int(255*x) for x in v[:3])

            columns_colors = []
            for col in df.columns:
                if col.startswith("smiles") or col.startswith("identifier") or col.startswith("molecule"):
                    columns_colors += ["#FFFFFF"]
                else:
                    columns_colors += [model_colors.get(col.split('_')[-1], "#FFFFFF")]
            column_colors = []

            for col_idx in range(num_cols):
                color = columns_colors[col_idx] if col_idx < len(columns_colors) else "#FFFFFF"
                col_format = workbook.add_format({'bg_color': color, 'align': 'left', 'valign': 'vcenter', 'bold': True, 'text_wrap': False})
                worksheet.write(0, col_idx, df.columns[col_idx], col_format)

            # Freeze the first row and apply autofilter
            worksheet.freeze_panes(1, 0)
            worksheet.autofilter(0, 0, len(df), num_cols - 1)

            # Add a second worksheet for the legend
            legend_ws = workbook.add_worksheet("Legend")
            legend_ws.write(0, 0, "model_id", workbook.add_format({'bold': True}))

            df = None
            for model_id in model_ids:
                slug = get_slug(model_id)
                dc = get_columns(model_id)
                data = {"model_id": [model_id]*len(dc), "slug": [slug]*len(dc), "col_num": [i+1 for i in range(len(dc))]}
                data = pd.DataFrame(data)
                data = pd.concat([data, dc], axis=1)
                if df is None:
                    df = data
                else:
                    df = pd.concat([df, data], axis=0).reset_index(drop=True)

            # Write DataFrame to legend worksheet
            for row_idx, row in df.iterrows():
                model_id = row["model_id"]
                color = model_colors.get(model_id, "#FFFFFF")
                cell_format = workbook.add_format({'bg_color': color, 'align': 'left', 'valign': 'vcenter'})
                legend_ws.write(row_idx + 1, 0, row["model_id"], cell_format)
                legend_ws.write(row_idx + 1, 1, row["slug"])
                legend_ws.write(row_idx + 1, 2, row["col_num"])
                legend_ws.write(row_idx + 1, 3, row["name"])
                legend_ws.write(row_idx + 1, 4, row["description"])
            # Set column headers
            legend_ws.write(0, 1, "slug", workbook.add_format({'bold': True}))
            legend_ws.write(0, 2, "col_num", workbook.add_format({'bold': True}))
            legend_ws.write(0, 3, "name", workbook.add_format({'bold': True}))
            legend_ws.write(0, 4, "description", workbook.add_format({'bold': True}))

            # Set column widths
            legend_ws.set_column('B:B', 20)  # slug column
            legend_ws.set_column('D:D', 30)  # name column
            legend_ws.set_column('E:E', 60)  # description column
            
            # Freeze first row
            legend_ws.freeze_panes(1, 0)

            writer.close()
            output.seek(0)

            st.download_button(
                label="📥 Download Excel file with molecules and predictions",
                data=output,
                file_name="merged_ersilia_file.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )