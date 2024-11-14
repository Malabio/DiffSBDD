#!/bin/bash

PROTEIN="LRRK2"
PROTEIN_VARIANT="G2019S"
PDB_ID="8TZC"
PDB_BINDING_POCKETS=("A:1885 A:1886 A:1887 A:1888 A:1889 A:1890 A:1893 A:1895 A:1902 A:1903 A:1904 A:1905 A:1906 A:1933 A:1947 A:1948 A:1949 A:1950 A:1951 A:1952 A:1953 A:1954 A:1955 A:1957 A:1994 A:1996 A:1998 A:1999 A:2000 A:2001 A:2016 A:2017 A:2019 A:2020 A:2033 A:2034 A:2035")
BASE_DIR="/home/markus/Malabio/DiffSBDD"

OPTIMIZATION_STEPS="150"
POPULATION_SIZE="8"
EVALUATION_STEPS="10"
TOP_K="2"
MOLS_TO_SAVE="ALL"

METRIC_CUTOFF=("-12 -10")
METRIC_SORTING=("QVinaScore$PROTEIN_VARIANT$PDB_ID QVinaScoreG2019S7LI3")
METRIC_ASCENDING=("False False")

PDBFILE="data/production/$PROTEIN/$PROTEIN_VARIANT/protein_only/$PDB_ID/$PDB_ID.pdb"
EVALUATION_FILE="results/$PROTEIN/evaluation/molecule_metrics.csv"
OUT_DIR="results/$PROTEIN/$PROTEIN_VARIANT/$PDB_ID/optimized$OPTIMIZATION_STEPS/"
CHECKPOINT="checkpoints/crossdocked_fullatom_cond.ckpt"

cd $BASE_DIR

#MAKE OUT DIR IF IT DOESN'T EXIST
mkdir -p $OUT_DIR

#Run the optimization file
echo "python optimize_results.py \
        --checkpoint $CHECKPOINT \
        --pdbfile $PDBFILE \
        --evaluation_file $EVALUATION_FILE \
        --outdir $OUT_DIR \
        --metric_sorting $METRIC_SORTING \
        --metric_cutoff $METRIC_CUTOFF \
        --metric_ascending $METRIC_ASCENDING \
        --variant $PROTEIN_VARIANT \
        --population_size $POPULATION_SIZE \
        --evolution_steps $EVALUATION_STEPS \
        --top_k $TOP_K \
        --timesteps $OPTIMIZATION_STEPS \
        --mols_to_save $MOLS_TO_SAVE \
        --resi_list $PDB_BINDING_POCKETS \
        --relax"
python optimize_results.py \
        --checkpoint $CHECKPOINT \
        --pdbfile $PDBFILE \
        --evaluation_file $EVALUATION_FILE \
        --outdir $OUT_DIR \
        --metric_sorting $METRIC_SORTING \
        --metric_cutoff $METRIC_CUTOFF \
        --metric_ascending $METRIC_ASCENDING \
        --variant $PROTEIN_VARIANT \
        --population_size $POPULATION_SIZE \
        --evolution_steps $EVALUATION_STEPS \
        --top_k $TOP_K \
        --timesteps $OPTIMIZATION_STEPS \
        --mols_to_save $MOLS_TO_SAVE \
        --resi_list $PDB_BINDING_POCKETS \
        --relax \


