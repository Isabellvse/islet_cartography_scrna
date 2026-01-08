#!/usr/bin/env bash
# ----------------------------- CONDA -----------------------------------
eval "$(conda shell.bash hook)"
conda activate /work/islet_cartography_scrna/scrna_cartography_gwas

# ----------------------------- USER CONFIG -----------------------------------
# Directories
WD="/work/islet_cartography_scrna"
GWAS_DIR="$WD/data/gwas"
ANN_DIR="$WD/data/anndata"
OUT_DIR="$GWAS_DIR/disease_scores"

# files
GS_FILE="$GWAS_DIR/gs_files/all_traits_geneset.gs"
COV_FILE="$GWAS_DIR/files/cov.tsv"
AD_OBJ="$ANN_DIR/AH_combined.h5ad"

# parameters
SPECIES="human"
ADJUST_PROP="manual_annotation"

# ----------------------------- DISEASE SCORES -----------------------------------
for trait in t2d_eur handedness_eur; do
    scdrs perform-downstream \
        --h5ad-file "$AD_OBJ" \
        --score-file "$OUT_DIR/${trait}.full_score.gz" \
        --out-folder "$OUT_DIR" \
        --group-analysis "$ADJUST_PROP" \
        --flag-filter-data True \
        --flag-raw-count False
done