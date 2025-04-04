import pandas as pd

def convert_xlsx_to_smi(xlsx_file, smiles_col, name_col=None, output_smi="list_of_mols.smi"):

    # Load Excel file
    df = pd.read_excel(xlsx_file)

    # Check if required columns exist
    if smiles_col not in df.columns:
        print(f"Error: Column '{smiles_col}' not found in the Excel file.")
        return

    # Extract required columns, keeping ligand_of_interest
    columns_to_keep = [smiles_col]
    if name_col and name_col in df.columns:
        columns_to_keep.append(name_col)

    df = df[columns_to_keep] # Keep relevant columns

    # Save to .smi format
    df[columns_to_keep].to_csv(output_smi, sep=" ", index=False, header=False)

# Example usage:
xlsx_file = "/Users/JB/Rotation_bkslab/250115_chaifold/20241209_mac1.xlsx"  # Replace with your actual file path
convert_xlsx_to_smi(xlsx_file, smiles_col="SMILES", name_col="Ligand_ID", output_smi="list_of_mols.smi")



def convert_col_to_csv(col_file, csv_file):

    # Read the .col file (Assuming it's space-separated or tab-separated)
    df = pd.read_csv(col_file, delim_whitespace=True, header=None)

    # Save as CSV
    df.to_csv(csv_file, index=False)

# Example usage:
col_file = "/Users/JB/Rotation_bkslab/250206_compare_ligands/max_tc_max_TC.col"  # Replace with your file path
csv_file = "ligand_check.csv"
convert_col_to_csv(col_file, csv_file)