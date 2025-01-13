from pymol import cmd

def extract_ligands(obj_name, cofactor_list):
    # Select all heteroatoms (which usually includes ligands)
    cmd.select("all_hetatm", f"{obj_name} and hetatm")

    # Create a selection string for excluding cofactors
    cofactor_selection = " or ".join([f"resn {cofactor}" for cofactor in cofactor_list])

    # Select ligands, excluding common cofactors
    cmd.select("ligands", f"all_hetatm and not ({cofactor_selection})")

    # Save the extracted ligands
    save_path = f"./{obj_name}_ligands.pdb"
    cmd.save(save_path, "ligands")

    # Optional: visualize the ligands
    cmd.show("sticks", "ligands")
    cmd.color("yellow", "ligands")
    cmd.zoom("ligands")

    # Clean up temporary selections
    cmd.delete("all_hetatm")
    cmd.delete("ligands")

# List of common cofactors to exclude
cofactor_list = ["HOH", "CL", "DMS", "TRS", "ATP", "ADP", "AMP", "GTP", "GDP", "FAD", "FMN", "NAD", "NADP", "COA", "HEME"]

# Get all objects in the current session
all_objects = cmd.get_object_list()

# Process each object
for obj in all_objects:
    extract_ligands(obj, cofactor_list)

print("Ligands extraction complete.")

