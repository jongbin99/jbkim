from pymol import cmd

# List of protein PDB files to load
protein_files = [
"7h7k.pdb",
"7h7j.pdb",
"7h7i.pdb",
"7h7g.pdb",
"7h7e.pdb",
"7h7c.pdb",
"7h7b.pdb",
"7h7a.pdb",
"7h78.pdb",
"7h77.pdb",
"7h76.pdb",
"7h75.pdb",
"7h73.pdb",
"7h72.pdb",
"7h71.pdb",
"7h6z.pdb",
"7h6y.pdb",
"7h6x.pdb",
"7h6w.pdb",
"7h6v.pdb",
"7h6q.pdb",
"7h6p.pdb",
"7h6n.pdb",
"7h6m.pdb",
"7h6l.pdb",
"7h6k.pdb",
"7h6j.pdb",
"7h7v.pdb",
"7h7s.pdb",
"7h7r.pdb",
"7h7o.pdb",
"7h7n.pdb",
"7h7l.pdb",
"7h8d.pdb",
"7h8c.pdb",
"7h8b.pdb",
"7h8a.pdb",
"7h89.pdb",
"7h87.pdb",
"7h86.pdb",
"7h85.pdb",
"7h84.pdb",
"7h83.pdb",
"7h82.pdb",
"7h81.pdb",
"7h80.pdb",
"7h7z.pdb",
"7h7y.pdb",
"7h7x.pdb",
"7h7w.pdb",
"7h8h.pdb",
"7h8f.pdb",
"7h95.pdb",
"7h94.pdb",
"7h93.pdb",
"7h92.pdb",
"7h91.pdb",
"7h90.pdb",
"7h8z.pdb",
"7h8y.pdb",
"7h8x.pdb",
"7h8w.pdb",
"7h8v.pdb",
"7h8u.pdb",
"7h8s.pdb",
"7h8r.pdb",
"7h8q.pdb",
"7h8p.pdb",
"7h8o.pdb",
"7h8n.pdb",
"7h8m.pdb",
"7h8k.pdb",
"7h8j.pdb",
"7h9j.pdb",
"7h9i.pdb",
"7h9h.pdb",
"7h9g.pdb",
"7h9f.pdb",
"7h9e.pdb",
"7h9d.pdb",
"7h9b.pdb",
"7h9a.pdb",
"7h99.pdb",
"7h98.pdb",
"7h97.pdb",
"7h96.pdb",
    # Add more protein files as needed
]

# List of known cofactors to exclude
cofactors = ['TRS', 'DMS', 'NAD', 'FAD', 'HEM', 'COA', 'FMN', 'ATP']

def is_ligand(resn):
    """Check if the residue is a ligand and not a cofactor."""
    return resn not in cofactors

def select_chains_with_ligands(protein):
    """Select chains that contain ligands."""
    chains_with_ligands = []
    cmd.select('ligands', f'{protein} and organic')
    ligands = cmd.get_model('ligands')
    
    for ligand in ligands.atom:
        if is_ligand(ligand.resn):
            chains_with_ligands.append(ligand.chain)
    
    return list(set(chains_with_ligands))

def remove_unaligned_chains(protein, aligned_chains):
    """Remove chains that are not aligned."""
    all_chains = cmd.get_chains(protein)
    for chain in all_chains:
        if chain not in aligned_chains:
            cmd.remove(f'{protein} and chain {chain}')

# Load all protein files
for protein_file in protein_files:
    cmd.load(protein_file)

# Perform alignment
reference_protein = protein_files[0].split('.')[0]  # Use the first protein as the reference
reference_chains = select_chains_with_ligands(reference_protein)

aligned_chains = {reference_protein: reference_chains}

for protein_file in protein_files[1:]:
    protein = protein_file.split('.')[0]
    chains = select_chains_with_ligands(protein)
    aligned_chains[protein] = []

    for ref_chain in reference_chains:
        for chain in chains:
            alignment_score = cmd.align(f'{protein} and chain {chain}', f'{reference_protein} and chain {ref_chain}')
            if alignment_score[0] < 2.0:  # Adjust the threshold as needed
                aligned_chains[protein].append(chain)
                break

# Remove unaligned chains
for protein, chains in aligned_chains.items():
    remove_unaligned_chains(protein, chains)

# Save the session
cmd.save('aligned_proteins.pse')

