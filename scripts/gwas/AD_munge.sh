#!/usr/bin/env bash
# ----------------------------- INFO -----------------------------------------


# ----------------------------- CONDA -----------------------------------------
eval "$(conda shell.bash hook)"
conda activate "/work/islet_cartography_scrna/scrna_cartography_gwas"

# ----------------------------- USER CONFIG -----------------------------------
# Define folders
WD="/work/islet_cartography_scrna"
GWAS_DIR="$WD/data/gwas"
PZ_DIR="$GWAS_DIR/pval_zscore"
GS_DIR="$GWAS_DIR/gs_files"

# ----------------------------- MUNGE -----------------------------------------
# Select top 1,00 genes and use z-score weights
scdrs munge-gs \
    --out-file "${GS_DIR}/all_traits_geneset.gs" \
    --zscore-file "${PZ_DIR}/zstat.tsv" \
    --fdr 0.05 \
    --n-max 100

# ----------------------------- CLEANUP ------------------------------------
RESULT_FILES=("${GS_DIR}/all_traits_geneset.gs")

ALL_EXIST=true
for f in "${RESULT_FILES[@]}"; do
    [ -f "$f" ] || { ALL_EXIST=false; break; }
done

if $ALL_EXIST; then
    rm -f *.tsv
    echo "Finished munging, you can find output file in ${GS_DIR}/all_traits_geneset.gs"
fi