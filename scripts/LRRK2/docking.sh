#!/bin/bash

PROTEIN="LRRK2"
PDB_IDS=("8TZC" "7LI3" "8TXZ")
PROTEIN_VARIANTS=("G2019S" "G2019S" "wild_type")
PDB_BINDING_POCKETS=("A:1885 A:1886 A:1887 A:1888 A:1889 A:1893 A:1895 A:1904 A:1905 A:1906 A:1933 A:1947 A:1948 A:1949 A:1950 A:1951 A:1952 A:1953 A:1954 A:1955 A:1957 A:1994 A:1996 A:1997 A:1998 A:1999 A:2000 A:2001 A:2002 A:2003 A:2014 A:2015 A:2016 A:2017 A:2019 A:2020"
                    "A:1885 A:1886 A:1887 A:1888 A:1889 A:1890 A:1893 A:1895 A:1902 A:1903 A:1904 A:1905 A:1906 A:1933 A:1947 A:1948 A:1949 A:1950 A:1951 A:1952 A:1953 A:1954 A:1955 A:1957 A:1994 A:1996 A:1998 A:1999 A:2000 A:2001 A:2016 A:2017 A:2019 A:2020 A:2033 A:2034 A:2035"
                    "A:1885 A:1886 A:1887 A:1888 A:1889 A:1890 A:1891 A:1892 A:1893 A:1895 A:1902 A:1903 A:1904 A:1905 A:1906 A:1933 A:1947 A:1948 A:1949 A:1950 A:1951 A:1952 A:1953 A:1955 A:1956 A:1957 A:1994 A:1996 A:1998 A:1999 A:2000 A:2001 A:2002 A:2003 A:2016 A:2017 A:2018 A:2019 A:2020 A:2033 A:2034 A:2035")
BASE_DIR="/home/markus/Malabio/DiffSBDD"
OPTIMIZATION_STEPS=("100")

cd $BASE_DIR

# Loop over the PDB IDs and queue the docking jobs.
counter2=0
while [ $counter2 -lt ${#PDB_IDS[@]} ]; do
    variant_sdf=${PROTEIN_VARIANTS[$counter2]}
    pdb_id_sdf=${PDB_IDS[$counter2]}
    counter=0
    pdb_length=${#PROTEIN_VARIANTS[@]}
    while [ $counter -lt $pdb_length ]; do
        # Defining the dirs
        pdb_variant=${PROTEIN_VARIANTS[$counter]}
        pdb_id=${PDB_IDS[$counter]}
        binding_pocket=${PDB_BINDING_POCKETS[$counter]}
        SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/generated/sdf/"
        OUT_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/generated/docked/$pdb_variant/$pdb_id/"
        PDBQT_DIR="$BASE_DIR/data/production/$PROTEIN/$pdb_variant/pdbqt/$pdb_id/"
        # Make out dir if it doesn't exist
        mkdir -p $OUT_DIR
        # echo the python command to be able to see the progress
        echo "python analysis/docking.py --pdbqt_dir $PDBQT_DIR --sdf_dir $SDF_DIR --out_dir $OUT_DIR --new_receptor_target $pdb_id --new_target_binding_pocket $binding_pocket"
        # Dock the ligands
        python analysis/docking.py \
                --pdbqt_dir $PDBQT_DIR \
                --sdf_dir $SDF_DIR \
                --out_dir $OUT_DIR \
                --new_receptor_target $pdb_id \
                --new_target_binding_pocket $binding_pocket 
        ((counter++))
    done
    # Also Dock the optimized molecules.
    counter=0
    while [ $counter -lt $pdb_length ]; do
        for opt_steps in "${OPTIMIZATION_STEPS[@]}"; do
            # Defining the dirs
            pdb_variant=${PROTEIN_VARIANTS[$counter]}
            pdb_id=${PDB_IDS[$counter]}
            binding_pocket=${PDB_BINDING_POCKETS[$counter]}
            SDF_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/optimized$opt_steps/sdf/"
            OUT_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/optimized$opt_steps/docked/$pdb_variant/$pdb_id/"
            PDBQT_DIR="$BASE_DIR/data/production/$PROTEIN/$pdb_variant/pdbqt/$pdb_id/"
            # Make out dir if it doesn't exist
            mkdir -p $OUT_DIR
            # echo the python command to be able to see the progress
            echo "python analysis/docking.py --pdbqt_dir $PDBQT_DIR --sdf_dir $SDF_DIR --out_dir $OUT_DIR --new_receptor_target $pdb_id --new_target_binding_pocket $binding_pocket"
            # Dock the ligands
            python analysis/docking.py \
                    --pdbqt_dir $PDBQT_DIR \
                    --sdf_dir $SDF_DIR \
                    --out_dir $OUT_DIR \
                    --new_receptor_target $pdb_id \
                    --new_target_binding_pocket $binding_pocket 
            ((counter++))
        done
    done
    # Also dock for the known inhibitors
    SDF_DIR="$BASE_DIR/results/$PROTEIN/inhibitors/"
    OUT_DIR="$BASE_DIR/results/$PROTEIN/$variant_sdf/$pdb_id_sdf/inhibitors/docked/"
    PDBQT_DIR="$BASE_DIR/data/production/$PROTEIN/$variant_sdf/pdbqt/$pdb_id_sdf/"
    pdb_variant=${PROTEIN_VARIANTS[$counter2]}
    pdb_id=${PDB_IDS[$counter2]}
    binding_pocket=${PDB_BINDING_POCKETS[$counter2]}
    #Make dir if it doesn't exist
    mkdir -p $OUT_DIR
    # echo the python command to be able to see the progress
    echo "python analysis/docking.py --pdbqt_dir $PDBQT_DIR --sdf_dir $SDF_DIR --out_dir $OUT_DIR --new_receptor_target $pdb_id_sdf --new_target_binding_pocket $binding_pocket"
    # Dock the ligands
    python analysis/docking.py \
            --pdbqt_dir $PDBQT_DIR \
            --sdf_dir $SDF_DIR \
            --out_dir $OUT_DIR \
            --new_receptor_target $pdb_id \
            --new_target_binding_pocket $binding_pocket 
    ((counter2++))
done