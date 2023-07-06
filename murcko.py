# Prototyping of the function to obtain the Standard Murco Scaffold from a list of SMILES.
# It is expected that the list is present, one molecule per line, in the file 'input.txt'.
# For testing purposes only. Eventually to be merged into a functioning tool.
#SkunkVersion 2023-07-06-1615.UTC

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# Read the input text file
with open('input.txt', 'r') as file:
    smiles_list = file.read().splitlines()

# Process each SMILES string
for smiles in smiles_list:
    # Convert SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)

    # Generate the Murcko scaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)

    # Convert the scaffold to a SMILES string
    scaffold_smiles = Chem.MolToSmiles(scaffold)

    # Print the Murcko scaffold SMILES
    print("SMILES:", smiles)
    print("Murcko Scaffold:", scaffold_smiles)
    print()
