#!/bin/bash

# Define the common variables
PROTEIN="LRRK2"
VARIANTS=("G2019S" "G2019S" "wild_type")
PDB_IDS=("8TZC" "7LI3" "8TXZ")
BASE_DIR="/home/markus/Malabio/DiffSBDD"
MODEL="crossdocked"

# Activate conda environment
conda init
conda activate mgltools

# Change to the base directory
cd $BASE_DIR/analysis

# Loop over the PDB IDs and create the pdbqt files from the pdb files
counter=0
while [ $counter -lt ${#PDB_IDS[@]} ]; do
    VARIANT=${VARIANTS[$counter]}
    PDB_ID=${PDB_IDS[$counter]}
    PDBDIR="$BASE_DIR/data/production/$PROTEIN/$VARIANT/protein_only/$PDB_ID/"
    PDBQT_DIR="$BASE_DIR/data/production/$PROTEIN/$VARIANT/pdbqt/$PDB_ID/"
    # Make out dir if it doesn't exist
    mkdir -p $PDBQT_DIR
    echo "Creating pdbqt files for $PROTEIN $VARIANT $PDB_ID..."
    echo "PDBDIR: $PDBDIR"
    echo "PDBQT_DIR: $PDBQT_DIR"
    python docking_py27.py $PDBDIR $PDBQT_DIR $MODEL
    ((counter++))
done