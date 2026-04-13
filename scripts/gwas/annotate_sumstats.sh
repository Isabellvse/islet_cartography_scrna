#!/usr/bin/env bash

# ----------------------------- USER CONFIG -----------------------------------
WD="/work/islet_cartography_scrna"
MAGMA_PROG="$WD/MAGMA"
REF_DIR="$WD/data/gwas/ref_magma"
GWAS_DIR="$WD/data/gwas"
SUMS_DIR="$GWAS_DIR/sumstats_processed"
OUT_DIR="$GWAS_DIR/sumstats_annotated"

cd "${MAGMA_PROG}"

# ----------------------------- LOOP ------------------------------------------
# Skip header (NR>1)
awk -F'\t' 'NR>1 {print $1,$2,$5,$6}' "$GWAS_DIR/files/magma_process_file.txt" | \
while read -r GENE_ANNOT GWAS_PVAL REF_BIM GENE_LOC; do

    BASENAME=$(basename "${GENE_ANNOT}" .genes.annot)
    OUT_PREFIX="${OUT_DIR}/${BASENAME}"

    [[ -f "${REF_DIR}/${REF_BIM}" ]] || { echo "Missing ${REF_BIM}"; exit 1; }
    [[ -f "${REF_DIR}/${GENE_LOC}" ]] || { echo "Missing ${GENE_LOC}"; exit 1; }
    [[ -f "${SUMS_DIR}/${GWAS_PVAL}" ]] || { echo "Missing ${GWAS_PVAL}"; exit 1; }

    echo "Processing ${BASENAME}"

    # --- Generate prefix for reference build ----
    REF_PREFIX="${REF_BIM%.bim}"
    
    # --- Annotation step ---
    ./magma --annotate window=10,10 \
        --snp-loc "${REF_DIR}/${REF_BIM}" \
        --gene-loc "${REF_DIR}/${GENE_LOC}" \
        --out "${OUT_PREFIX}"

    # --- Gene analysis step ---
    
    ./magma --bfile "${REF_DIR}/${REF_PREFIX}" \
        --pval "${SUMS_DIR}/${GWAS_PVAL}" use='SNP,P' ncol='N' \
        --gene-annot "${OUT_PREFIX}.genes.annot" \
        --out "${OUT_PREFIX}"

done
