import argparse
from pathlib import Path

import torch
import ast
from rdkit import Chem
import numpy as np

import pandas as pd

from analysis.metrics import MoleculeProperties
from utils import read_sdf_file

def make_metrics_df(mol_dir, target):
    # Load the docking results
    docking_result_path = Path(mol_dir, f'qvina_scores.csv')
    docking_result = pd.read_csv(docking_result_path, index_col=0)
    if len(docking_result) == 0:
        print(f"No docking results found for {target} and molecule directory {mol_dir}")
        return None
    # Get the first metrics for the molecules
    filenames = docking_result["ligand"]
    scores = docking_result["scores"]
    generated = [True if not 'knowninhib' in filename else False for filename in filenames]
    target_variant = [target for filename in filenames]
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
    lipinski = []
    lipinski_violations = []
    # Define the mol metrics calculator
    mol_metrics = MoleculeProperties()
    # Loop through the molecules
    for i, file in enumerate(filenames):
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
            lipinski.append(None)
            lipinski_violations.append([None, None, None, None, None, None])
            continue
        # Get the metrics for the molecules
        smiles.append(Chem.MolToSmiles(mol))
        qed.append(mol_metrics.calculate_qed(mol))
        sa.append(mol_metrics.calculate_sa(mol))
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
        f'QVinaScore{target}': scores,
        'QED': qed,
        'SA': sa,
        'LogP': logp,
        'Weight': weight,
        'HDonor': hdonor,
        'HAcceptor': hacceptor,
        'RotatableBonds': rotatable,
        'TPSA': tpsa,
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

def merge_metric_dfs(all_dfs, targets):
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

    
    # final_df = None
    # for i, df in enumerate(all_dfs):
    #     if df is None:
    #         continue
    #     rec_name = targets[i]
    #     if final_df is None:
    #         final_df = df
    #     else:
    #         # Get common columns
    #         # Check if the filenames are already in the dataframe
    #         df_filenames = df['Filename'].values
    #         metrics_df_filenames = final_df['Filename'].values
    #         common_filenames = list(set(df_filenames).intersection(set(metrics_df_filenames)))
    #         if len(common_filenames) == 0:
    #             if np.array([c in df.columns for c in final_df.columns]).all():
    #                 final_df = pd.concat([final_df, df], ignore_index=True)
    #             else:
    #                 # If no common filenames, merge on all common columns
    #                 common_columns = list(set(final_df.columns).intersection(set(df.columns)))
    #                 # Merge the dataframes based on common columns other than 'Smile' and 'QVinaScore'
    #                 merged_df = pd.merge(final_df, df, on=common_columns, how='outer')
    #                 final_df = merged_df
    #         else:
    #             # If the filenames are already in the dataframe, find the column where values are different, this will be the QVinaScore column that needs to be added.
    #             # Check if the values are not the same for the columns
    #             # Get the correct values from the metrics_df for the filenames
    #             metrics_df_common = final_df.loc[final_df["Filename"].isin(common_filenames)]
    #             # Make sure they are in the same order by sorting on the filenames and idx
    #             metrics_df_common = metrics_df_common.sort_values(by=["Filename"])
    #             df = df.sort_values(by=["Filename"])
    #             final_df = final_df.sort_values(by=["Filename"])
    #             for col in df.columns:
    #                 # Check if the column is not in the metrics_df because then this is a new QVinaScore column.
    #                 if col not in final_df.columns:
    #                     final_df[col] = None
    #                     final_df.loc[(final_df["Filename"].isin(common_filenames)), col] = df[col]                    # Convert to string to be able to compare NaN values as well
    #                 elif not np.all(metrics_df_common[col].values.astype(str) == df[col].values.astype(str)):
    #                     if col == "Smile":
    #                         continue
    #                     elif col in ["QED", "SA", "LogP", "Weight", "HDonor", "HAcceptor", "RotatableBonds", "TPSA", "Lipinski"]:
    #                         # Take Average
    #                         final_df.loc[final_df["Filename"].isin(common_filenames), col] = (final_df.loc[final_df["Filename"].isin(common_filenames), col] + df[col]) / 2
    #                     elif "QVina" in col:
    #                         # If the values are not the same we know it is because they are NaN in the metrics_df as the QVinaScore is not present in the metrics_df yet
    #                         final_df.loc[final_df["Filename"].isin(common_filenames), col] = df[col]
    #                     else:
    #                         final_df.loc[final_df["Filename"].isin(common_filenames), col] = df[col]
    #                 else:
    #                     pass
    # return final_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mol_dirs', type=str, nargs='+', default=None)
    parser.add_argument('--docking_targets', type=str, nargs='+', default=None)
    parser.add_argument('--outdir', type=Path, default=None)
    args = parser.parse_args()

    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # Getting mol_dirs docking targets
    docking_targets = []
    for d in args.mol_dirs:
        d = d.split("/")
        for i in range(len(d)):
            if d[-1-i] in args.docking_targets:
                docking_targets.append(d[-1-i])
                break
    all_dfs = []
    for j, mol_dir in enumerate(args.mol_dirs):
        # Get the metrics for all molecule directories
        result_df = make_metrics_df(mol_dir, docking_targets[j])
        # Append the dataframe to the list
        all_dfs.append(result_df)

    # Merge the dataframes
    final_df = merge_metric_dfs(all_dfs, docking_targets)

    #Order by docking score
    final_df = final_df.sort_values(by=[f'QVinaScore'+pdb for pdb in np.unique(docking_targets)], ascending=True)

    # Save the dataframe
    final_df.to_csv(args.outdir / 'molecule_metrics.csv', index=False)




    #     sdf_files = list(Path(mol_dir).glob('[!.]*.sdf'))
    #     # Add inhibitor sdf files
    #     receptor_name_gen = mol_dir.split("/")
    #     results_idx = [i for i, s in enumerate(receptor_name_gen) if 'results' in s][0]
    #     receptor_name_gen = receptor_name_gen[results_idx+2]
    #     for i, docking_target in enumerate(args.docking_targets):
    #         docking_result_path = Path(mol_dir, f'qvina_scores.csv')
    #         docking_result = pd.read_csv(docking_result_path)
    #         if len(docking_result) == 0:
    #             print(f"No docking results found for {docking_target} and molecule directory {mol_dir}")
    #             continue
    #         filenames = docking_result.iloc[:,2]
    #         generated = [True if not 'knowninhib' in filename else False for filename in filenames]
    #         docking_result["generated"] = generated
    #         filenames = []
    #         scores = []
    #         generated = []
    #         for i, row in docking_result.iterrows():
    #             score = row["scores"]
    #             filename = row["ligand"]
    #             gen = row["generated"]
    #             score = score.replace("nan", "None")
    #             score = ast.literal_eval(score)
    #             scores.append(score)
    #             filenames.append(filename)
    #             generated.append(gen)

    #         new_docking_results = []
    #         sdf_indexes = []
    #         new_generated = []
    #         new_filenames = []

    #         # Find valid molecules
    #         total_molecules = 0
    #         valid_molecules = []
    #         for k, filename in enumerate(filenames):
    #             molecules = read_sdf_file(Path(filename), sanitize=True)
    #             # This happens when a file is read where multiple conformations of the same molecule are present.
    #             if len(molecules) > len(scores[k]):
    #                 assert len(scores[k]) == 1, "Multiple conformations of the same molecule are present in the file, but multiple scores are present as well."
    #                 molecules = [molecules[0]]
    #             total_molecules += len(molecules)
    #             score = scores[k]
    #             gen = generated[k]
    #             if gen == False and len(score) == 1:
    #                 score = [score[0] for i in range(len(molecules))]
    #             #Do not include directory in filename
    #             name = filename.split("/")[-1]
    #             for i, mol in enumerate(molecules):
    #                 if mol is not None:
    #                     valid_molecules.append(mol)
    #                     sdf_indexes.append(i)
    #                     new_generated.append(gen)
    #                     new_filenames.append(filename)
    #                     new_docking_results.append(score[i])

    #         smiles = [Chem.MolToSmiles(mol) for mol in valid_molecules]
        
    #         print(f"Found {len(valid_molecules)}/{total_molecules} valid molecules")

    #         mol_metrics = MoleculeProperties()

    #         valid_molecules = [valid_molecules]

    #         all_qed, all_sa, all_logp, all_lipinski, all_lipinski_violations, all_lipinski_values = mol_metrics.evaluate(valid_molecules, return_lipinski_violations=True, return_diversity=False)
    #         # all_qed, all_sa, all_logp, all_lipinski, per_pocket_diversity = mol_metrics.evaluate(valid_molecules, return_lipinski_violations=True, return_diversity=False)

    #         #Getting values out of the lipinski vlaues.
    #         all_weight = [x[0] for x in all_lipinski_values[0]]
    #         all_hdonor = [x[1] for x in all_lipinski_values[0]]
    #         all_hacceptor = [x[2] for x in all_lipinski_values[0]]
    #         all_rotatable = [x[4] for x in all_lipinski_values[0]]
    #         all_tpsa = [x[5] for x in all_lipinski_values[0]]

    #         #Save metrics as csv
    #         metrics = { f'Smile': smiles, 
    #                     f'QVinaScore{docking_target}': new_docking_results, 
    #                     'sdfFileIdx':sdf_indexes,
    #                     'Filename': new_filenames,
    #                     'Generated': new_generated,
    #                     'QED': all_qed[0], 
    #                     'SA': all_sa[0], 
    #                     'LogP': all_logp[0],
    #                     'Weight': all_weight,
    #                     'HDonor': all_hdonor,
    #                     'HAcceptor': all_hacceptor,
    #                     'RotatableBonds': all_rotatable,
    #                     'Lipinski': all_lipinski[0],
    #                     'TPSA': all_tpsa, 
    #                     'Lipinski_violation_rule1': [x[0] for x in all_lipinski_violations[0]],
    #                     'Lipinski_violation_rule2': [x[1] for x in all_lipinski_violations[0]],
    #                     'Lipinski_violation_rule3': [x[2] for x in all_lipinski_violations[0]],
    #                     'Lipinski_violation_rule4': [x[3] for x in all_lipinski_violations[0]],
    #                     'Lipinski_violation_rule5': [x[4] for x in all_lipinski_violations[0]],
    #                     'Lipinski_violation_rule6': [x[5] for x in all_lipinski_violations[0]],
    #                     #    'Diversity': per_pocket_diversity[0],
    #                     }
    #         metrics_df = pd.DataFrame(metrics)
    #         metrics_df["ReceptorGen"] = receptor_name_gen
    #         all_metric_dfs.append(metrics_df)
    #         all_receptor_names_gen.append(receptor_name_gen)
    #         all_receptor_names_dock.append(docking_target)

    #     # Do it for the known inhibitors
    #     docking_target = receptor_name_gen
    #     # Get idx of first receptor_name_gen in mol_dir
    #     print(mol_dir)
    #     print(docking_target)
    #     inhib_dir = mol_dir.split("/")
    #     variant_idx = [i for i, s in enumerate(inhib_dir) if docking_target in s][0]
    #     inhib_dir = "/".join(inhib_dir[:variant_idx+2])
    #     docking_result_path = Path(inhib_dir, f'inhibitors/qvina2_scores.csv')
    #     docking_result = pd.read_csv(docking_result_path)
    #     filenames = docking_result.iloc[:,2]
    #     generated = [True if not 'knowninhib' in filename else False for filename in filenames]
    #     docking_result["generated"] = generated
    #     filenames = []
    #     scores = []
    #     generated = []
    #     for i, row in docking_result.iterrows():
    #         score = row["scores"]
    #         filename = row["ligand"]
    #         gen = row["generated"]
    #         score = score.replace("nan", "None")
    #         score = ast.literal_eval(score)
    #         scores.append(score)
    #         filenames.append(filename)
    #         generated.append(gen)

    #     new_docking_results = []
    #     sdf_indexes = []
    #     new_generated = []
    #     new_filenames = []

    #     # Find valid molecules
    #     valid_molecules = []
    #     for k, filename in enumerate(filenames):
    #         molecules = read_sdf_file(Path(filename), sanitize=True)
    #         score = scores[k]
    #         gen = generated[k]
    #         if gen == False and len(score) == 1:
    #             score = [score[0] for i in range(len(molecules))]
    #         #Do not include directory in filename
    #         name = filename.split("/")[-1]
    #         for i, mol in enumerate(molecules):
    #             if mol is not None:
    #                 valid_molecules.append(mol)
    #                 sdf_indexes.append(i)
    #                 new_generated.append(gen)
    #                 new_filenames.append(filename)
    #                 new_docking_results.append(score[i])

    #     smiles = [Chem.MolToSmiles(mol) for mol in valid_molecules]

    #     print(f"Found {len(valid_molecules)}/{len(molecules)} valid molecules")

    #     mol_metrics = MoleculeProperties()

    #     valid_molecules = [valid_molecules]

    #     all_qed, all_sa, all_logp, all_lipinski, all_lipinski_violations, all_lipinski_values = mol_metrics.evaluate(valid_molecules, return_lipinski_violations=True, return_diversity=False)
    #     # all_qed, all_sa, all_logp, all_lipinski, per_pocket_diversity = mol_metrics.evaluate(valid_molecules, return_lipinski_violations=True, return_diversity=False)

    #     #Getting values out of the lipinski vlaues.
    #     all_weight = [x[0] for x in all_lipinski_values[0]]
    #     all_hdonor = [x[1] for x in all_lipinski_values[0]]
    #     all_hacceptor = [x[2] for x in all_lipinski_values[0]]
    #     all_rotatable = [x[4] for x in all_lipinski_values[0]]
    #     all_tpsa = [x[5] for x in all_lipinski_values[0]]

    #     #Save metrics as csv
    #     metrics = { f'Smile': smiles, 
    #                 f'QVinaScore{docking_target}': new_docking_results, 
    #                 'sdfFileIdx':sdf_indexes,
    #                 'Filename': new_filenames,
    #                 'Generated': new_generated,
    #                 'QED': all_qed[0], 
    #                 'SA': all_sa[0], 
    #                 'LogP': all_logp[0],
    #                 'Weight': all_weight,
    #                 'HDonor': all_hdonor,
    #                 'HAcceptor': all_hacceptor,
    #                 'RotatableBonds': all_rotatable,
    #                 'TPSA': all_tpsa,
    #                 'Lipinski': all_lipinski[0], 
    #                 'Lipinski_violation_rule1': [x[0] for x in all_lipinski_violations[0]],
    #                 'Lipinski_violation_rule2': [x[1] for x in all_lipinski_violations[0]],
    #                 'Lipinski_violation_rule3': [x[2] for x in all_lipinski_violations[0]],
    #                 'Lipinski_violation_rule4': [x[3] for x in all_lipinski_violations[0]],
    #                 'Lipinski_violation_rule5': [x[4] for x in all_lipinski_violations[0]],
    #                 'Lipinski_violation_rule6': [x[5] for x in all_lipinski_violations[0]],
    #                 #    'Diversity': per_pocket_diversity[0],
    #                 }
    #     metrics_df = pd.DataFrame(metrics)
    #     metrics_df["ReceptorGen"] = receptor_name_gen
    #     all_metric_dfs.append(metrics_df)
    #     all_receptor_names_gen.append(receptor_name_gen)
    #     all_receptor_names_dock.append(docking_target)

    # # First merge the dataframes which have molecules generated for the same receptor
    # metrics_df = {f'{rec_name}': None for rec_name in np.unique(all_receptor_names_gen)}
    # for i, df in enumerate(all_metric_dfs):
    #     rec_name = all_receptor_names_gen[i]
    #     if metrics_df[rec_name] is None:
    #         metrics_df[rec_name] = df
    #     else:
    #         # Get common columns
    #         # Check if the filenames are already in the dataframe
    #         df_filenames = df['Filename'].values
    #         metrics_df_filenames = metrics_df[rec_name]['Filename'].values
    #         common_filenames = list(set(df_filenames).intersection(set(metrics_df_filenames)))
    #         if len(common_filenames) == 0:
    #             # If no common filenames, merge on all common columns
    #             common_columns = list(set(metrics_df[rec_name].columns).intersection(set(df.columns)))
    #             # Merge the dataframes based on common columns other than 'Smile' and 'QVinaScore'
    #             merged_df = pd.merge(metrics_df[rec_name], df, on=common_columns, how='outer')
    #             metrics_df[rec_name] = merged_df
    #         else:
    #             # If the filenames are already in the dataframe, find the column where values are different, this will be the QVinaScore column that needs to be added.
    #             # Check if the values are not the same for the columns
    #             # Get the correct values from the metrics_df for the filenames
    #             metrics_df_common = metrics_df[rec_name].loc[metrics_df[rec_name]["Filename"].isin(common_filenames)]
    #             # Make sure they are in the same order by sorting on the filenames and idx
    #             metrics_df_common = metrics_df_common.sort_values(by=["Filename", "sdfFileIdx"])
    #             df = df.sort_values(by=["Filename", "sdfFileIdx"])
    #             metrics_df[rec_name] = metrics_df[rec_name].sort_values(by=["Filename", "sdfFileIdx"])
    #             for col in df.columns:
    #                 # Check if the column is not in the metrics_df because then this is a new QVinaScore column.
    #                 if col not in metrics_df[rec_name].columns:
    #                     metrics_df[rec_name][col] = df[col]
    #                 # Convert to string to be able to compare NaN values as well
    #                 elif not np.all(metrics_df_common[col].values.astype(str) == df[col].values.astype(str)):
    #                     # Sort the values based on the filenames
    #                     metrics_df[rec_name] = metrics_df[rec_name].sort_values(by="Filename")
    #                     # If the values are not the same we know it is because they are NaN in the metrics_df as the QVinaScore is not present in the metrics_df yet
    #                     metrics_df[rec_name].loc[metrics_df[rec_name]["Filename"].isin(common_filenames), col] = df[col]
    #                 else:
    #                     pass


    #         # merge_columns = list(set(metrics_df[rec_name].columns).intersection(set(df.columns)))
    #         # merge_columns = [col for col in merge_columns if "QVinaScore" not in col]
    #         # # Merge the dataframes based on common columns other than 'Smile' and 'QVinaScore'
    #         # merged_df = pd.merge(metrics_df[rec_name], df, on=merge_columns, how='outer')
    #         # metrics_df[rec_name] = merged_df