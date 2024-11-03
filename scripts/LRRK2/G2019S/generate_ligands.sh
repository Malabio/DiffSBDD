#!/bin/bash

# Define the common variables
VARIANT="G2019S"
PDB_ID="7LI3"
PROTEIN="LRRK2"
BASE_DIR="/home/markus/Malabio/DiffSBDD"

CKPT="checkpoints/crossdocked_fullatom_cond.ckpt"
PDBFILE="data/production/$PROTEIN/$VARIANT/protein_only/$PDB_ID.pdb"
OUTDIR="results/$PROTEIN/$VARIANT/$PDB_ID/generated/"
RESI_LIST="A:1885 A:1886 A:1887 A:1888 A:1889 A:1890 A:1891 A:1892 A:1893 A:1895 A:1902 A:1903 A:1904 A:1905 A:1906 A:1933 A:1947 A:1948 A:1949 A:1950 A:1951 A:1952 A:1953 A:1955 A:1956 A:1957 A:1994 A:1996 A:1998 A:1999 A:2000 A:2001 A:2002 A:2003 A:2016 A:2017 A:2018 A:2019 A:2020 A:2033 A:2034 A:2035"
BATCH_SIZE=48
N_SAMPLES=48
COMMON_OPTS="--pdbfile $PDBFILE --resi_list $RESI_LIST --batch_size $BATCH_SIZE --n_samples $N_SAMPLES --sanitize --relax --all_frags"
MIN_NODES=25
MAX_NODES=70

# Change to the base directory
cd $BASE_DIR

# Create the output directory if it does not exist recursively
mkdir -p $OUTDIR

# Range of num_nodes_lig values (from 25 to 70)
for NUM_NODES in $(seq $MIN_NODES $MAX_NODES); do
    OUTFILE="$OUTDIR/${VARIANT}_${PDB_ID}_size${NUM_NODES}.sdf"
    echo "Running ligand generation with num_nodes_lig=$NUM_NODES..."
    python generate_ligands.py $CKPT $COMMON_OPTS --outfile $OUTFILE --num_nodes_lig $NUM_NODES
done

echo "Batch ligand generation complete."
