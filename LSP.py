# Prototyping of the function to identify the longest shortest path (LSP) of a molecule
# For testing purposes only. Eventually to be merged into a functioning tool which combines the LSP
# and the Beemis-Murcko scaffold (MBM).
# SkunkVersion 2023-07-06-1615.UTC

import sys
from rdkit import Chem
import numpy as np

# Check if the input file is provided as a command line argument
if len(sys.argv) < 2:
    print("Please provide the input file as a command line argument.")
    sys.exit(1)

# Read the input file
input_file = sys.argv[1]
smiles_list = []

with open(input_file, 'r') as file:
    smiles_list = file.read().splitlines()

# Process each SMILES string
for smiles in smiles_list:
    # Convert SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)

    # Get the number of atoms in the molecule
    num_atoms = mol.GetNumAtoms()

    # Create an adjacency matrix for the molecule
    adj_matrix = Chem.GetAdjacencyMatrix(mol)

    # Initialize the distance matrix with large values
    dist_matrix = np.full((num_atoms, num_atoms), np.inf)

    # Set the diagonal elements to 0
    np.fill_diagonal(dist_matrix, 0)

    # Initialize the path matrix to store the atom indices
    path_matrix = np.zeros((num_atoms, num_atoms), dtype=object)

    # Update the distance matrix with the initial connectivity information
    for i in range(num_atoms):
        for j in range(num_atoms):
            if adj_matrix[i, j] == 1:
                dist_matrix[i, j] = 1
                path_matrix[i, j] = [i, j]

    # Use the Floyd-Warshall algorithm to calculate the shortest paths
    for k in range(num_atoms):
        for i in range(num_atoms):
            for j in range(num_atoms):
                if dist_matrix[i, k] + dist_matrix[k, j] < dist_matrix[i, j]:
                    dist_matrix[i, j] = dist_matrix[i, k] + dist_matrix[k, j]
                    path_matrix[i, j] = path_matrix[i, k] + path_matrix[k, j][1:]

    # Find the longest shortest path in the molecule
    longest_path_index = np.unravel_index(np.argmax(dist_matrix), dist_matrix.shape)
    longest_path = path_matrix[longest_path_index]

    # Generate the SMILES string for the longest shortest path
    longest_path_atoms = longest_path
    longest_path_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=longest_path_atoms)

    print("Original SMILES:", smiles)
    print("Longest shortest path:", longest_path)
    print("Longest shortest path length:", len(longest_path) - 1)
    print("Longest shortest path SMILES:", longest_path_smiles)
    print()
