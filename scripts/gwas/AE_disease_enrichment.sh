#!/usr/bin/env bash
# ----------------------------- CONDA -----------------------------------------
eval "$(conda shell.bash hook)"
conda activate /work/islet_cartography_scrna/scrna_cartography_gwas

# ----------------------------- USER CONFIG -----------------------------------
# Directories
WD="/work/islet_cartography_scrna"
GWAS_DIR="$WD/data/gwas"
ANN_DIR="$WD/data/anndata"
OUT_DIR="$GWAS_DIR/disease_scores"

# files
GS_FILE="$GWAS_DIR/gs_files/pops_genes.gs"
COV_FILE="$GWAS_DIR/files/cov.tsv"
AD_OBJ="$ANN_DIR/AH_combined.h5ad"

# parameters
SPECIES="human"
ADJUST_PROP="manual_annotation"

# ----------------------------- DISEASE SCORES -----------------------------------
scdrs compute-score \
    --h5ad-file "$AD_OBJ" \
    --h5ad-species "$SPECIES" \
    --gs-file "$GS_FILE" \
    --gs-species "$SPECIES" \
    --cov-file "$COV_FILE" \
    --flag-filter-data False \
    --flag-raw-count False \
    --flag-return-ctrl-raw-score True \
    --flag-return-ctrl-norm-score True \
    --out-folder "$OUT_DIR"