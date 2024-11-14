#!/bin/bash

PROTEIN="LRRK2"
PDB_IDS=("8TZC" "7LI3" "8TXZ")
PROTEIN_VARIANTS=("G2019S" "G2019S" "wild_type")
BASE_DIR="/home/markus/Malabio/DiffSBDD"

OPTIMIZATION_STEPS=("100" "50")

# Go to the base directory
cd $BASE_DIR

counter=0
while [ $counter -lt ${#PDB_IDS[@]} ]; do
    variant_sdf=${PROTEIN_VARIANTS[$counter]}
    pdb_id_sdf=${PDB_IDS[$counter]}
    counter2=0
    while [ $counter2 -lt ${#PDB_IDS[@]} ]; do
        variant=${PROTEIN_VARIANTS[$counter2]}
        pdb_id=${PDB_IDS[$counter2]}
        SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/generated/docked/$variant/$pdb_id/"
        # Make the Qvina score files
        python make_qvina_scores_files.py \
            --sdf_dir $SDF_DIR \
            --receptor_file $BASE_DIR/data/production/$PROTEIN/$variant/pdbqt/$pdb_id/$pdb_id.pdbqt 
        for opt_steps in "${OPTIMIZATION_STEPS[@]}"; do
            SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/optimized$opt_steps/docked/$variant/$pdb_id/"
            # Make the Qvina score files
            python make_qvina_scores_files.py \
                --sdf_dir $SDF_DIR \
                --receptor_file $BASE_DIR/data/production/$PROTEIN/$variant/pdbqt/$pdb_id/$pdb_id.pdbqt
        done
        ((counter2++))
    done
    SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/inhibitors/docked/"
    # Make the Qvina score files
    python make_qvina_scores_files.py \
        --sdf_dir $SDF_DIR \
        --receptor_file $BASE_DIR/data/production/$PROTEIN/$variant_sdf/pdbqt/$pdb_id_sdf/$pdb_id_sdf.pdbqt
    ((counter++))
done