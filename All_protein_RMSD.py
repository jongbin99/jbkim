import os
import csv
from pymol import cmd

# Set directory containing PDB files
input_dir = "/Users/JB/Rotation_bkslab/250115_chaifold/Organized_CIFS/idx_0_pdb_final/all_files"
output_csv = "rmsd_results.csv"

if not os.path.exists(input_dir):
    print(f"[ERROR] Directory not found: {input_dir}")
    exit()

reference_files = []
predicted_files = []
pocket_residues = "21+22+23+48+49+52+125+126+129+130+154+155+156+157+160" #residues in the binding pocket

# Get reference and predicted files
for file_name in os.listdir(input_dir):
    if file_name.endswith(".pdb"):
        if "_pred_chainA.pdb" in file_name:
            predicted_files.append(file_name)
        elif file_name.startswith("mac-x"):
            reference_files.append(file_name)

rmsd_results = []

# Iterate over reference files
for ref_file in reference_files:
    prefix = os.path.splitext(ref_file)[0]
    pred_file = next((f for f in predicted_files if f.startswith(prefix) and "_pred_chainA" in f), None)

    if pred_file is None:
        print(f"[WARNING] No matching predicted file for {ref_file}. Skipping.")
        continue

    ref_path = os.path.join(input_dir, ref_file)
    pred_path = os.path.join(input_dir, pred_file)

    if not os.path.exists(ref_path) or not os.path.exists(pred_path):
        print(f"[ERROR] Missing file: {ref_path} or {pred_path}")
        continue

    ref_object = f"{prefix}_R"
    pred_object = f"{prefix}_P"
    
    try:
        cmd.load(ref_path, ref_object)
        cmd.load(pred_path, pred_object)
        print(f"[INFO] Loaded {ref_object} and {pred_object}")
    except Exception as e:
        print(f"[ERROR] Failed to load structures: {e}")
        continue

    ref_count = cmd.count_atoms(ref_object)
    pred_count = cmd.count_atoms(pred_object)

    if ref_count == 0 or pred_count == 0:
        print(f"[ERROR] No atoms found in {ref_object} or {pred_object}. Skipping.")
        continue

#Align and compute total RMSD
    try:
        alignment_result = cmd.align(pred_object, ref_object)
        full_rmsd = alignment_result[0]
        print(f"[INFO] Full Structure RMSD for {prefix}: {full_rmsd:.3f} Å")
    except Exception as e:
        print(f"[ERROR] Alignment failed: {e}")
        continue

#Compute pocket RMSD
    try:
        pocket_ref = f"{ref_object} and resi {pocket_residues}"
        pocket_pred = f"{pred_object} and resi {pocket_residues}"
        pocket_rmsd = cmd.rms_cur(pocket_pred, pocket_ref)
        print(f"[INFO] Pocket Side-Chain RMSD for {prefix}: {pocket_rmsd:.3f} Å")

    except Exception as e:
        print(f"[ERROR] Failed to compute pocket RMSD for {prefix}: {e}")
        pocket_rmsd = "N/A"
    
#Compute pocket side-chain RMSD only
    try:
        sc_ref = f"{ref_object} and resi {pocket_residues} and not name CA+C+N+O"
        sc_pred = f"{pred_object} and resi {pocket_residues} and not name CA+C+N+O"
        side_chain_rmsd = cmd.rms_cur(sc_pred, sc_ref)
        print(f"[INFO] Pocket Side-Chain RMSD for {prefix}: {side_chain_rmsd:.3f} Å")

    except Exception as e:
        print(f"[ERROR] Failed to compute pocket RMSD for {prefix}: {e}")
        side_chain_rmsd = "N/A"
        
#Compute pocket backbone RMSD only
    try:
        backbone_selection_ref = f"{ref_object} and resi {pocket_residues} and name CA+C+N+O"
        backbone_selection_pred = f"{pred_object} and resi {pocket_residues} and name CA+C+N+O"
        ref_count_bb = cmd.count_atoms(backbone_selection_ref)
        pred_count_bb = cmd.count_atoms(backbone_selection_pred)

        if ref_count_bb == 0 or pred_count_bb == 0:
            print(f"[WARNING] No backbone atoms in {prefix}. Skipping.")
            backbone_rmsd = "N/A"
        else:
            backbone_rmsd = cmd.rms_cur(backbone_selection_pred, backbone_selection_ref)
            print(f"[INFO] Backbone RMSD for {prefix}: {backbone_rmsd:.3f} Å")
    except Exception as e:
        print(f"[ERROR] Backbone RMSD failed: {e}")
        backbone_rmsd = "N/A"

    rmsd_results.append([prefix, full_rmsd, pocket_rmsd, side_chain_rmsd, backbone_rmsd])
    cmd.delete(ref_object)
    cmd.delete(pred_object)

# Save results to CSV
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Protein", "Full Structure RMSD (Å)", "Pocket RMSD (Å)", "Side-Chain RMSD (Å)", "Backbone RMSD (Å)"])
    writer.writerows(rmsd_results)

print(f"[INFO] RMSD results saved to {output_csv}")