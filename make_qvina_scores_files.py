import argparse
import os
import pandas as pd
import rdkit.Chem as Chem
import re
from pathlib import Path

def make_qvina_scores_file(sdf_dir, receptor_file):
    # Get all sdf files in the directory
    sdf_files = [f for f in os.listdir(sdf_dir) if f.endswith('.sdf')]

    # Create a list to store the scores
    scores = []
    files = []

    # Loop through each sdf file
    for sdf_file in sdf_files:
        # Read the sdf file
        # print(f'Reading {sdf_file}...')
        suppl = Chem.SDMolSupplier(os.path.join(sdf_dir, sdf_file))
        qvina = 10000
        # Get all the QVina score Remark field
        for mol in suppl:
            if mol is not None:
                if mol.HasProp('REMARK'):
                    remark = mol.GetProp('REMARK')
                    # Get the Vina score
                    match = re.search(r"VINA RESULT:\s+(-?\d+\.\d+)", remark)
                    # Print only the extracted number
                    vina_score = float(match.group(1)) if match else None
                    if vina_score is not None and vina_score < qvina:
                        qvina = vina_score
        # Append the min score to the list
        scores.append(qvina)
        files.append(Path(sdf_dir, sdf_file))
    receptors = [receptor_file] * len(sdf_files)
    # Make the DataFrame
    df = pd.DataFrame({'receptor': receptors, 'ligand': files, 'scores': scores})
    # Save the DataFrame to a CSV file
    df.to_csv(Path(sdf_dir, 'qvina_scores.csv'), index=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a CSV file with QVina scores for a directory of sdf files.')
    parser.add_argument('--sdf_dirs', type=str, nargs='+', default=None, help='The directory containing the sdf files.')
    parser.add_argument('--receptor_file', type=str, default='ligand', help='The ligand name to use in the output file.')
    args = parser.parse_args()

    for sdf_dir in args.sdf_dirs:
        # Make sdf_dir if it does not exist
        if not os.path.exists(sdf_dir):
            os.makedirs(sdf_dir)
        make_qvina_scores_file(sdf_dir, args.receptor_file)