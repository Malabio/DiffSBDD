#!/bin/bash

PROTEIN="LRRK2"
DOCKING_TARGETS=("G2019S" "wild_type")
BASE_DIR="/home/markus/Malabio/DiffSBDD"
OUTDIR="$BASE_DIR/results/$PROTEIN/evaluation/"
MOL_DIRS=()
OPTIMIZATION_STEPS=("100")

for target1 in "${DOCKING_TARGETS[@]}"; do
    for target2 in "${DOCKING_TARGETS[@]}"; do
        MOL_DIRS+=("$BASE_DIR/results/$PROTEIN/$target1/docked/$target2/generated/")
        for opt_steps in "${OPTIMIZATION_STEPS[@]}"; do
            MOL_DIRS+=("$BASE_DIR/results/$PROTEIN/$target1/docked/$target2/optimized$opt_steps/")
        done
    done
done

# Go to the base directory
cd $BASE_DIR

# Make output directory if it doesn't exist
mkdir -p $OUTDIR

echo "Running evaluation for $PROTEIN"
echo "Output directory: $OUTDIR"
echo "Molecule directories: ${MOL_DIRS[@]}"

Run the evaluate script
python evaluate.py \
    --outdir $OUTDIR \
    --mol_dirs "${MOL_DIRS[@]}" \
    --docking_targets "${DOCKING_TARGETS[@]}"