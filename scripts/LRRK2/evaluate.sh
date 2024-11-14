#!/bin/bash

PROTEIN="LRRK2"
PDB_IDS=("8TZC" "7LI3" "8TXZ")
PROTEIN_VARIANTS=("G2019S" "G2019S" "wild_type")
BASE_DIR="/home/markus/Malabio/DiffSBDD"
OUTDIR="$BASE_DIR/results/$PROTEIN/evaluation/"
MOL_DIRS=()
SDF_VARIANTS=()
SDF_IDS=()
DOCKING_VARIANTS=()
DOCKING_IDS=()
OPTIMIZATION_STEPS=("100" "50")

counter=0
while [ $counter -lt ${#PDB_IDS[@]} ]; do
    variant_sdf=${PROTEIN_VARIANTS[$counter]}
    pdb_id_sdf=${PDB_IDS[$counter]}
    counter2=0
    while [ $counter2 -lt ${#PDB_IDS[@]} ]; do
        variant=${PROTEIN_VARIANTS[$counter2]}
        pdb_id=${PDB_IDS[$counter2]}
        MOL_DIRS+=("$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/generated/docked/$variant/$pdb_id/")
        SDF_VARIANTS+=("$variant_sdf")
        SDF_IDS+=("$pdb_id_sdf")
        DOCKING_VARIANTS+=("$variant")
        DOCKING_IDS+=("$pdb_id")
        for opt_steps in "${OPTIMIZATION_STEPS[@]}"; do
            MOL_DIRS+=("$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/optimized$opt_steps/docked/$variant/$pdb_id/")
            SDF_VARIANTS+=("$variant_sdf")
            SDF_IDS+=("$pdb_id_sdf")
            DOCKING_VARIANTS+=("$variant")
            DOCKING_IDS+=("$pdb_id")
        done
        ((counter2++))
    done
    MOL_DIRS+=("$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/inhibitors/docked/")
    SDF_VARIANTS+=("inhibitors")
    SDF_IDS+=("inhibitors")
    DOCKING_VARIANTS+=("$variant_sdf")
    DOCKING_IDS+=("$pdb_id_sdf")
    ((counter++))
done

# Go to the base directory
cd $BASE_DIR

# Make output directory if it doesn't exist
mkdir -p $OUTDIR

echo "Running evaluation for $PROTEIN"
echo "Output directory: $OUTDIR"
echo "Molecule directories: ${MOL_DIRS[@]}"
echo "SDF variants: ${SDF_VARIANTS[@]}"
echo "SDF IDs: ${SDF_IDS[@]}"
echo "Docking variants: ${DOCKING_VARIANTS[@]}"
echo "Docking IDs: ${DOCKING_IDS[@]}"

Run the evaluate script
python evaluate.py \
    --outdir $OUTDIR \
    --mol_dirs "${MOL_DIRS[@]}" \
    --sdf_variants "${SDF_VARIANTS[@]}" \
    --sdf_pdb_ids "${SDF_IDS[@]}" \
    --docking_variants "${DOCKING_VARIANTS[@]}" \
    --docking_pdb_ids "${DOCKING_IDS[@]}"