# ToBeFixed: The counting of similars above a certain threshold is wrong - it appears to be counting itself. So that needs to fixed.
# Additional features to add: separating out exacts for non-redundant output (optional)
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity

def load_molecules_from_sdf(sdf_file):
    supplier = Chem.SDMolSupplier(sdf_file)
    molecules = [mol for mol in supplier if mol is not None]
    return molecules

def calculate_morgan_fingerprints(molecules, radius=2, nBits=2048):
    fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits) for mol in molecules]
    return fingerprints

def calculate_tanimoto_matrix(fingerprints):
    n = len(fingerprints)
    tanimoto_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            tanimoto_matrix[i, j] = TanimotoSimilarity(fingerprints[i], fingerprints[j])
            tanimoto_matrix[j, i] = tanimoto_matrix[i, j]
    return tanimoto_matrix

def find_nearest_neighbors(tanimoto_matrix):
    n = len(tanimoto_matrix)
    closest_neighbors = [np.argsort(-row)[1] for row in tanimoto_matrix]
    closest_neighbors_sim = [tanimoto_matrix[i, closest_neighbors[i]] for i in range(n)]
    return closest_neighbors, closest_neighbors_sim

def count_neighbors_above_thresholds(tanimoto_matrix, thresholds):
    counts = {threshold: [] for threshold in thresholds}
    for row in tanimoto_matrix:
        for threshold in thresholds:
            counts[threshold].append(np.sum(row > threshold))
    return counts

def main(sdf_file, csv_output, sdf_output, matrix_output):
    molecules = load_molecules_from_sdf(sdf_file)
    fingerprints = calculate_morgan_fingerprints(molecules)
    tanimoto_matrix = calculate_tanimoto_matrix(fingerprints)

    smiles_list = [Chem.MolToSmiles(mol) for mol in molecules]
    names_list = [mol.GetProp('_Name') if mol.HasProp('_Name') else '' for mol in molecules]

    closest_neighbors, closest_neighbors_sim = find_nearest_neighbors(tanimoto_matrix)
    closest_neighbors_smiles = [smiles_list[idx] for idx in closest_neighbors]
    closest_neighbors_names = [names_list[idx] for idx in closest_neighbors]

    thresholds = [0.5, 0.7, 0.8, 0.9]
    neighbors_count = count_neighbors_above_thresholds(tanimoto_matrix, thresholds)

    # Write to CSV (summary)
    data = {
        'Name': names_list,
        'SMILES': smiles_list,
        'Closest Neighbor Name': closest_neighbors_names,
        'Closest Neighbor SMILES': closest_neighbors_smiles,
        'Closest Neighbor Similarity': closest_neighbors_sim,
        'Neighbors >0.5': neighbors_count[0.5],
        'Neighbors >0.7': neighbors_count[0.7],
        'Neighbors >0.8': neighbors_count[0.8],
        'Neighbors >0.9': neighbors_count[0.9]
    }
    df = pd.DataFrame(data)
    df.to_csv(csv_output, index=False)

    # Write the complete Tanimoto similarity matrix to CSV
    tanimoto_df = pd.DataFrame(tanimoto_matrix, index=names_list, columns=names_list)
    tanimoto_df.to_csv(matrix_output)

    # Write to SDF
    writer = Chem.SDWriter(sdf_output)
    for i, mol in enumerate(molecules):
        mol.SetProp("Closest Neighbor Name", closest_neighbors_names[i])
        mol.SetProp("Closest Neighbor SMILES", closest_neighbors_smiles[i])
        mol.SetProp("Closest Neighbor Similarity", str(closest_neighbors_sim[i]))
        mol.SetProp("Neighbors >0.5", str(neighbors_count[0.5][i]))
        mol.SetProp("Neighbors >0.7", str(neighbors_count[0.7][i]))
        mol.SetProp("Neighbors >0.8", str(neighbors_count[0.8][i]))
        mol.SetProp("Neighbors >0.9", str(neighbors_count[0.9][i]))
        writer.write(mol)
    writer.close()

# Example usage
sdf_file = 'input.sdf'  # Replace with your input SDF file path
csv_output = 'output_summary.csv'
sdf_output = 'output.sdf'
matrix_output = 'tanimoto_matrix.csv'
main(sdf_file, csv_output, sdf_output, matrix_output)
