#!/usr/bin/env bash
# ----------------------------- CONDA -----------------------------------
eval "$(conda shell.bash hook)"
conda activate /work/islet_cartography_scrna/scrna_cartography_gwas

# ----------------------------- USER CONFIG -----------------------------------
# Directories
WD="/work/islet_cartography_scrna"
GWAS_DIR="$WD/data/gwas"
ANN_DIR="$WD/data/anndata"
SCORE_DIR="$GWAS_DIR/disease_scores"
OUT_DIR="$GWAS_DIR/downstream_analysis"

# files
GS_FILE="$GWAS_DIR/gs_files/all_traits_geneset.gs"
COV_FILE="$GWAS_DIR/files/cov.tsv"
AD_OBJ="$ANN_DIR/AH_combined.h5ad"

# parameters
SPECIES="human"
GROUP="manual_annotation"

# ----------------------------- DISEASE SCORES -----------------------------------
# Get files files
FILES="$SCORE_DIR/*full_score.gz"

for z in $FILES; do
    scdrs perform-downstream \
        --h5ad-file "$AD_OBJ" \
        --score-file "$z" \
        --out-folder "$OUT_DIR" \
        --group-analysis "$GROUP" \
        --flag-filter-data True \
        --flag-raw-count False
done