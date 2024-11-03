#!/bin/bash

PROTEIN="LRRK2"
PROTEIN_VARIANTS=("wild_type" "G2019S")
PDB_BINDING_POCKETS=("A:1885 A:1886 A:1887 A:1888 A:1889 A:1893 A:1895 A:1904 A:1905 A:1906 A:1933 A:1947 A:1948 A:1949 A:1950 A:1951 A:1952 A:1953 A:1954 A:1955 A:1957 A:1994 A:1996 A:1997 A:1998 A:1999 A:2000 A:2001 A:2002 A:2003 A:2014 A:2015 A:2016 A:2017 A:2019 A:2020"
                    "A:1885 A:1886 A:1887 A:1888 A:1889 A:1890 A:1893 A:1895 A:1902 A:1903 A:1904 A:1905 A:1906 A:1933 A:1947 A:1948 A:1949 A:1950 A:1951 A:1952 A:1953 A:1954 A:1955 A:1957 A:1994 A:1996 A:1998 A:1999 A:2000 A:2001 A:2016 A:2017 A:2019 A:2020 A:2033 A:2034 A:2035"
                    )
BASE_DIR="/home/markus/Malabio/DiffSBDD"
OPTIMIZATION_STEPS=("100")

cd $BASE_DIR

# Loop over the PDB IDs and queue the docking jobs.
counter2=0
for variant in "${PROTEIN_VARIANTS[@]}"; do
    counter=0
    pdb_length=${#PROTEIN_VARIANTS[@]}
    while [ $counter -lt $pdb_length ]; do
        # Defining the dirs
        pdb_id_pdbfile=${PROTEIN_VARIANTS[$counter]}
        binding_pocket=${PDB_BINDING_POCKETS[$counter]}
        SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant/generated/"
        OUT_DIR="$BASE_DIR/results/$PROTEIN/$variant/docked/$pdb_id_pdbfile/generated/"
        PDBQT_DIR="$BASE_DIR/data/production/$PROTEIN/$pdb_id_pdbfile/pdbqt/"
        # Make out dir if it doesn't exist
        mkdir -p $OUT_DIR
        # Dock the ligands
        python analysis/docking.py \
                --pdbqt_dir $PDBQT_DIR \
                --sdf_dir $SDF_DIR \
                --out_dir $OUT_DIR \
                --new_receptor_target $pdb_id_pdbfile \
                --new_target_binding_pocket $binding_pocket \
                --write_csv \
                --write_dict
        ((counter++))
    done
    # Also Dock the optimized molecules.
    counter=0
    while [ $counter -lt $pdb_length ]; do
        for opt_steps in "${OPTIMIZATION_STEPS[@]}"; do
            # Defining the dirs
            pdb_id_pdbfile=${PROTEIN_VARIANTS[$counter]}
            binding_pocket=${PDB_BINDING_POCKETS[$counter]}
            SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant/optimized/"
            OUT_DIR="$BASE_DIR/results/$PROTEIN/$variant/docked/$pdb_id_pdbfile/optimized$opt_steps/"
            PDBQT_DIR="$BASE_DIR/data/production/$PROTEIN/$pdb_id_pdbfile/pdbqt/"
            # Make out dir if it doesn't exist
            mkdir -p $OUT_DIR
            # Dock the ligands
            python analysis/docking.py \
                    --pdbqt_dir $PDBQT_DIR \
                    --sdf_dir $SDF_DIR \
                    --out_dir $OUT_DIR \
                    --new_receptor_target $pdb_id_pdbfile \
                    --new_target_binding_pocket $binding_pocket \
                    --write_csv \
                    --write_dict
            ((counter++))
        done
    done
    # Also dock for the known inhibitors
    binding_pocket=${PDB_BINDING_POCKETS[$counter2]}
    SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant/inhibitors/"
    OUT_DIR="$BASE_DIR/results/$PROTEIN/$variant/docked/inhibitors/"
    PDBQT_DIR="$BASE_DIR/data/production/$PROTEIN/$variant/pdbqt/"
    pdb_id=${PROTEIN_VARIANTS[$counter2]}
    #Make dir if it doesn't exist
    mkdir -p $OUT_DIR
    python analysis/docking.py \
            --pdbqt_dir $PDBQT_DIR \
            --sdf_dir $SDF_DIR \
            --out_dir $OUT_DIR \
            --new_receptor_target $pdb_id \
            --new_target_binding_pocket $binding_pocket \
            --write_csv \
            --write_dict
    ((counter2++))
done