import os
import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.drawOptions.addAtomIndices = True

excel_file = "/Users/JB/Rotation_bkslab/250203_alphafold3/20241209_mac1.xlsx" #The file with list of SMILES
ref_dir = "/Users/JB/Rotation_bkslab/250203_alphafold3/PDBS_lig" #Directory with experimental structures
pred_dir = "/Users/JB/Rotation_bkslab/Mac1_docked_poses/aligned_ligands_dock" #Directory with predicted structures
output_csv = "/Users/JB/Rotation_bkslab/Mac1_docked_poses/docked_L-RMSD.csv" #Output CSV file

df = pd.read_excel(excel_file)
df['Dataset_ID'] = df['Dataset_ID'].astype(str)
df['SMILES'] = df['SMILES'].astype(str)

rmsd_results = []

# Loop through reference files
for ref_file in os.listdir(ref_dir):
    if not ref_file.endswith(".pdb"):
        continue

    # Extract the prefix
    ref_prefix = os.path.splitext(ref_file)[0]

    # Find the corresponding predicted file
    pred_file = f"{ref_prefix}.pdb"
    pred_path = os.path.join(pred_dir, pred_file)
    ref_path = os.path.join(ref_dir, ref_file)

    if not os.path.exists(pred_path):
        print(f"[WARNING] Predicted file not found for {ref_prefix}. Skipping.")
        continue

    # Retrieve SMILES from the Excel file
    tmplt_smiles_row = df.loc[df['Dataset_ID'] == ref_prefix, 'SMILES']
    if tmplt_smiles_row.empty or pd.isna(tmplt_smiles_row.values[0]):
        print(f"[WARNING] No SMILES found for {ref_prefix}. Skipping.")
        continue
    tmplt_smiles = str(tmplt_smiles_row.values[0])

    # Convert the SMILES to an RDKit molecule
    tmplt_mol = Chem.MolFromSmiles(tmplt_smiles)
    if not tmplt_mol:
        print(f"[ERROR] Invalid SMILES for {ref_prefix}: {tmplt_smiles}. Skipping.")
        continue

    # Load reference and predicted structures
    reference = Chem.MolFromPDBFile(ref_path, removeHs=False)
    predicted = Chem.MolFromPDBFile(pred_path, removeHs=False)

    if not reference or not predicted:
        print(f"[ERROR] Failed to load structures for {ref_prefix}. Skipping.")
        continue

    # Assign bond orders
    try:
        reference = AllChem.AssignBondOrdersFromTemplate(tmplt_mol, reference)
        predicted = AllChem.AssignBondOrdersFromTemplate(tmplt_mol, predicted)
    except Exception as e:
        print(f"[ERROR] Bond order assignment failed for {ref_prefix}: {e}")
        continue

    # Perform substructure matching and renumber atoms
    matches = predicted.GetSubstructMatches(reference)
    if not matches:
        print(f"[ERROR] No substructure matches found for {ref_prefix}. Skipping.")
        continue
    match_list = matches[0]
    predicted = Chem.RenumberAtoms(predicted, match_list)

    # Display molecular structures
    Draw.MolsToGridImage([reference, predicted], legends=["Reference", "Predicted"], molsPerRow=2)

    # Find RMSD
    try:
        rmsd = AllChem.GetBestRMS(reference, predicted)
        print(f"RMSD for {ref_prefix}: {rmsd:.3f} Å")
        rmsd_results.append([ref_prefix, rmsd])
    except Exception as e:
        print(f"[ERROR] RMSD calculation failed for {ref_prefix}: {e}")
        continue

# Write RMSD results to a CSV file
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["ID", "RMSD"])
    writer.writerows(rmsd_results)

print(f"[INFO] RMSD results saved to {output_csv}")