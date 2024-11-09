import argparse
from pathlib import Path
import openbabel as ob

import torch
from tqdm import tqdm
from rdkit import Chem
import numpy as np
import time

import pandas as pd

from analysis.metrics import MoleculeProperties
from utils import read_sdf_file

def make_metrics_df(mol_dir, sdf_variant, sdf_pdb_id, docking_variant, docking_pdb_id):
    # Load the docking results
    docking_result_path = Path(mol_dir, f'qvina_scores.csv')
    docking_result = pd.read_csv(docking_result_path, index_col=0)
    if len(docking_result) == 0:
        print(f"No docking results found for sdf_variant {sdf_variant} sdf_pdb_id {sdf_pdb_id}, docking_variant {docking_variant} docking_pdb_id {docking_pdb_id}, mol_dir {mol_dir}")
        return None
    
    # Get the first metrics for the molecules
    filenames = docking_result["ligand"]
    scores = docking_result["scores"]
    generated = [True if not 'knowninhib' in filename else False for filename in filenames]

    # Defining lists to save the rest
    smiles = []
    qed = []
    sa = []
    logp = []
    weight = []
    hdonor = []
    hacceptor = []
    rotatable = []
    tpsa = []
    pka_acid = []
    pka_base = []
    lipinski = []
    lipinski_violations = []

    # Define the mol metrics calculator
    mol_metrics = MoleculeProperties()

    # Loop through the molecules
    for i, file in tqdm(enumerate(filenames), total=len(filenames)):
        # Read the sdf file
        molecules = read_sdf_file(Path(file), sanitize=True)
        # As all molecules are the same in the files, just different conformation, we only need the first molecule
        mol = molecules[0]
        if mol is None:
            smiles.append(None)
            qed.append(None)
            sa.append(None)
            logp.append(None)
            weight.append(None)
            hdonor.append(None)
            hacceptor.append(None)
            rotatable.append(None)
            tpsa.append(None)
            pka_acid.append(None)
            pka_base.append(None)
            lipinski.append(None)
            lipinski_violations.append([None, None, None, None, None, None])
            continue
        
        # Get the metrics for the molecules
        smiles.append(Chem.MolToSmiles(mol))
        qed.append(mol_metrics.calculate_qed(mol))
        sa.append(mol_metrics.calculate_sa(mol))
        pka_acid_, pka_base_ = mol_metrics.calculate_pka(mol)
        pka_acid.append(pka_acid_)
        pka_base.append(pka_base_)
        lipinski_, lipinski_violations_, [weight_, hdonor_, hacceptor_, logp_, rotatable_, tpsa_] = mol_metrics.calculate_lipinski(mol, return_rule_violations=True)
        lipinski.append(lipinski_)
        lipinski_violations.append(lipinski_violations_)
        weight.append(weight_)
        hdonor.append(hdonor_)
        hacceptor.append(hacceptor_)
        logp.append(logp_)
        rotatable.append(rotatable_)
        tpsa.append(tpsa_)

    # Only have the last of the filenames path
    filenames = [file.split("/")[-1] for file in filenames]
    filenames = np.array(filenames)
    for i, f in enumerate(filenames):
        if "knowninhib" in f:
            f = f.split("_")[1:]
            filenames[i] = "target_" + "_".join(f)

    # Make Dataframe
    metrics = {
        'Smile': smiles,
        'Filename': filenames,
        'Generated': generated,
        f'QVinaScore{docking_variant}{docking_pdb_id}': scores,
        'QED': qed,
        'SA': sa,
        'LogP': logp,
        'Weight': weight,
        'HDonor': hdonor,
        'HAcceptor': hacceptor,
        'RotatableBonds': rotatable,
        'TPSA': tpsa,
        'pKaAcid': pka_acid,
        'pKaBase': pka_base,
        'Lipinski': lipinski,
        'Lipinski_violation_rule1': [x[0] for x in lipinski_violations],
        'Lipinski_violation_rule2': [x[1] for x in lipinski_violations],
        'Lipinski_violation_rule3': [x[2] for x in lipinski_violations],
        'Lipinski_violation_rule4': [x[3] for x in lipinski_violations],
        'Lipinski_violation_rule5': [x[4] for x in lipinski_violations],
        'Lipinski_violation_rule6': [x[5] for x in lipinski_violations],
    }
    metrics_df = pd.DataFrame(metrics)
    return metrics_df

def merge_metric_dfs(all_dfs):
    # Get all filenames and columns from the dataframes
    filenames = []
    columns = []
    for i, df in enumerate(all_dfs):
        if df is None:
            continue
        filenames += list(df['Filename'].values)
        columns += list(df.columns)
    # Get unique filenames and columns
    filenames = np.unique(filenames)
    columns = np.unique(columns)
    # Create the final dataframe
    final_df = pd.DataFrame(columns=columns)
    # Add the filenames
    final_df['Filename'] = filenames
    # Add a column counting number of times the filename is present
    final_df['FilenameCount'] = 0
    # Loop through the dataframes
    for i, df in enumerate(all_dfs):
        if df is None:
            continue
        for j, file in df.iterrows():
            filename_idx = np.where(final_df['Filename'] == file['Filename'])[0][0]
            filename_count = final_df.loc[final_df['Filename'] == file['Filename'], 'FilenameCount'].item()
            if filename_count == 0:
                final_df.loc[filename_idx, 'FilenameCount'] += 1
                final_df.loc[filename_idx, df.columns] = file
            else:
                final_df.loc[filename_idx, 'FilenameCount'] += 1
                for col in df.columns:
                    if col in ["QED", "SA", "LogP", "Weight", "HDonor", "HAcceptor", "RotatableBonds", "TPSA", "Lipinski"]:
                        # Take Average
                        final_df.loc[filename_idx, col] = final_df.loc[filename_idx, col] * filename_count / (filename_count + 1) + file[col] / (filename_count + 1)
                    if "QVinaScore" in col:
                        final_df.loc[filename_idx, col] = file[col]
                    else:
                        final_df.loc[filename_idx, col] = file[col]
    return final_df

def calculate_pka(df, protein):
    """
    This function calculates the pKa of the molecules in the dataframe and adds it to the dataframe.
    """

    # Define the mol metrics calculator
    mol_metrics = MoleculeProperties()

    # Define time, as only 500 requests per hour are allowed
    t1 = time.time()
    count = 0

    # Loop through the molecules
    for i, row in tqdm(df.iterrows(), total=len(df)):
        # Get the Filename
        filename = row['Filename']
        variant = filename.split("_")[0]
        pdb_id = filename.split("_")[1]
        # Get sdf name
        sdf_name = f"results/{protein}/{variant}/{pdb_id}/{row['Method']}/docked/{variant}/{pdb_id}/{filename}"
        # Get the molecule
        molecules = read_sdf_file(Path(sdf_name), sanitize=True)
        # As all molecules are the same in the files, just different conformation, we only need the first molecule
        mol = molecules[0]
        if mol is None:
            pka_acid = None
            pka_base = None
        else:
            # Get the metrics for the molecules
            t2 = time.time()
            if count == 0:
                count += 1
                pka_acid, pka_base = mol_metrics.calculate_pka(mol)
                continue
            if ((t2 - t1) / count) < (3600 / 500):
                # Sleep for the remaining time
                time.sleep(count*(3600 / 500) - (t2-t1))
            
            # Get the metrics for the molecules
            pka_acid, pka_base = mol_metrics.calculate_pka(mol)
            count += 1

        df.loc[i, "pKaAcid"] = pka_acid
        df.loc[i, "pKaBase"] = pka_base

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--evaluate_file', type=str, default=None)
    parser.add_argument('--protein', type=str, default=None)
    args = parser.parse_args()

    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # Load the evaluate file
    evaluate_file = pd.read_csv(args.evaluate_file, sep=',')

    # Merge the dataframes
    evaluate_pka = calculate_pka(evaluate_file, args.protein)

    # Save the dataframe
    evaluate_pka.to_csv(args.evaluate_file, index=False)