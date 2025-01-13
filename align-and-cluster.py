from pymol import cmd, stored
import numpy as np
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer

# List of protein PDB files to load
protein_files = [
    "7h7k.pdb", "7h7j.pdb", "7h7i.pdb", "7h7g.pdb", "7h7e.pdb", "7h7c.pdb", "7h7b.pdb", "7h7a.pdb",
    "7h78.pdb", "7h77.pdb", "7h76.pdb", "7h75.pdb", "7h73.pdb", "7h72.pdb", "7h71.pdb", "7h6z.pdb",
    "7h6y.pdb", "7h6x.pdb", "7h6w.pdb", "7h6v.pdb", "7h6q.pdb", "7h6p.pdb", "7h6n.pdb", "7h6m.pdb",
    "7h6l.pdb", "7h6k.pdb", "7h6j.pdb", "7h7v.pdb", "7h7s.pdb", "7h7r.pdb", "7h7o.pdb", "7h7n.pdb",
    "7h7l.pdb", "7h8d.pdb", "7h8c.pdb", "7h8b.pdb", "7h8a.pdb", "7h89.pdb", "7h87.pdb", "7h86.pdb",
    "7h85.pdb", "7h84.pdb", "7h83.pdb", "7h82.pdb", "7h81.pdb", "7h80.pdb", "7h7z.pdb", "7h7y.pdb",
    "7h7x.pdb", "7h7w.pdb", "7h8h.pdb", "7h8f.pdb", "7h95.pdb", "7h94.pdb", "7h93.pdb", "7h92.pdb",
    "7h91.pdb", "7h90.pdb", "7h8z.pdb", "7h8y.pdb", "7h8x.pdb", "7h8w.pdb", "7h8v.pdb", "7h8u.pdb",
    "7h8s.pdb", "7h8r.pdb", "7h8q.pdb", "7h8p.pdb", "7h8o.pdb", "7h8n.pdb", "7h8m.pdb", "7h8k.pdb",
    "7h8j.pdb", "7h9j.pdb", "7h9i.pdb", "7h9h.pdb", "7h9g.pdb", "7h9f.pdb", "7h9e.pdb", "7h9d.pdb",
    "7h9b.pdb", "7h9a.pdb", "7h99.pdb", "7h98.pdb", "7h97.pdb", "7h96.pdb",
]

# List of known cofactors to exclude
cofactors = ['NAD', 'FAD', 'HEM', 'COA', 'FMN', 'ATP', 'TRS', 'DMS']

# Define the number of clusters
n_clusters = 10  # Set your desired number of clusters here

def is_ligand(resn):
    """Check if the residue is a ligand and not a cofactor."""
    return resn not in cofactors

def get_ligand_positions(protein):
    """Get the positions of ligands in the protein, excluding cofactors."""
    cmd.select('ligands', f'{protein} and organic')
    ligands = cmd.get_model('ligands')
    
    positions = []
    for atom in ligands.atom:
        if is_ligand(atom.resn):
            positions.append([atom.coord[0], atom.coord[1], atom.coord[2]])
    
    return positions

def select_chains_with_ligands(protein):
    """Select chains that contain ligands, excluding cofactors."""
    chains_with_ligands = set()
    cmd.select('ligands', f'{protein} and organic')
    ligands = cmd.get_model('ligands')
    
    for ligand in ligands.atom:
        if is_ligand(ligand.resn):
            chains_with_ligands.add(ligand.chain)
    
    return list(chains_with_ligands)

def remove_unaligned_chains(protein, aligned_chains):
    """Remove chains that are not aligned."""
    all_chains = cmd.get_chains(protein)
    for chain in all_chains:
        if chain not in aligned_chains:
            cmd.remove(f'{protein} and chain {chain}')

def pad_positions(positions, max_length):
    """Pad positions with NaNs to ensure they all have the same length."""
    padded = np.full((max_length, 3), np.nan)
    padded[:len(positions), :] = positions
    return padded.flatten()

# Load all protein files
for protein_file in protein_files:
    cmd.load(protein_file)

# Identify chains with ligands and align based on those chains
reference_protein = protein_files[0].split('.')[0]  # Use the first protein as the reference
reference_chains = select_chains_with_ligands(reference_protein)

aligned_chains = {reference_protein: reference_chains}

for protein_file in protein_files[1:]:
    protein = protein_file.split('.')[0]
    chains = select_chains_with_ligands(protein)
    aligned_chains[protein] = []
    
    for ref_chain in reference_chains:
        for chain in chains:
            alignment_score = cmd.align(f'{protein} and chain {chain} and not solvent', f'{reference_protein} and chain {ref_chain} and not solvent')
            if alignment_score[0] < 2.0:  # Adjust the threshold as needed
                aligned_chains[protein].append(chain)
                break

# Remove unaligned chains but keep waters
for protein, chains in aligned_chains.items():
    remove_unaligned_chains(protein, chains)

# Extract ligand positions from the aligned chains, excluding cofactors
all_positions = []
chain_ligand_positions = {}
max_length = 0
for protein_file in protein_files:
    protein = protein_file.split('.')[0]
    positions = []
    chain_positions = {}
    for chain in aligned_chains[protein]:
        chain_pos = get_ligand_positions(f'{protein} and chain {chain}')
        chain_positions[chain] = chain_pos
        positions.extend(chain_pos)
    all_positions.append(positions)
    chain_ligand_positions[protein] = chain_positions
    max_length = max(max_length, len(positions))

# Pad all positions to the same length
all_positions_padded = [pad_positions(positions, max_length) for positions in all_positions]

# Convert to numpy array
all_positions_padded = np.array(all_positions_padded)

# Impute missing values with 0
imputer = SimpleImputer(missing_values=np.nan, strategy='constant', fill_value=0)
all_positions_imputed = imputer.fit_transform(all_positions_padded)

# Perform clustering with the user-defined number of clusters
kmeans = KMeans(n_clusters=n_clusters)
labels = kmeans.fit_predict(all_positions_imputed)

# Find the centroids of each cluster
centroids = kmeans.cluster_centers_

# Assign chains with multiple ligands based on the closest ligand to the cluster centroids
for i, protein_file in enumerate(protein_files):
    protein = protein_file.split('.')[0]
    chain_positions = chain_ligand_positions[protein]
    for chain, positions in chain_positions.items():
        if len(positions) > 1:  # Chain contains more than one ligand
            min_dist = float('inf')
            closest_cluster = labels[i]
            for pos in positions:
                pos_array = np.array(pos).reshape(1, -1)
                for cluster_id, centroid in enumerate(centroids):
                    dist = np.linalg.norm(pos_array - centroid[:3])
                    if dist < min_dist:
                        min_dist = dist
                        closest_cluster = cluster_id
            labels[i] = closest_cluster

# Create and save separate sessions for each cluster
for cluster_id in range(n_clusters):
    # Select all proteins in the current cluster
    cmd.select(f'cluster_{cluster_id}', 'none')
    for i, protein_file in enumerate(protein_files):
        if labels[i] == cluster_id:
            protein = protein_file.split('.')[0]
            cmd.select(f'cluster_{cluster_id}', f'cluster_{cluster_id} or {protein}')
    # Save the cluster to a separate session file
    cmd.save(f'cluster_{cluster_id}.pse', f'cluster_{cluster_id}')

print(f"Clustering complete with {n_clusters} clusters. Each cluster has been saved to a separate session file.")

