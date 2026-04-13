#!/usr/bin/env bash

# ----------------------------- USER CONFIG -----------------------------------
WD="/work/islet_cartography_scrna"
REF_DIR="$WD/data/gwas/ref_magma"
GWAS_DIR="$WD/data/gwas"
SUMS_DIR="$GWAS_DIR/sumstats_annotated"
OUT_DIR="$GWAS_DIR/gs_files"

# ----------------------------- LOOP ------------------------------------------
# Output file
mkdir -p "$OUT_DIR"
echo -e "TRAIT\tGENESET" > "${OUT_DIR}/all_traits_geneset.gs"

awk -F'\t' 'NR>1 {print $1,$6}' "$GWAS_DIR/files/magma_process_file.txt" | \
while read -r GENE_ANNOT GENE_LOC; do

    BASENAME=$(basename "${GENE_ANNOT}" .genes.annot)
    
    echo "Running ${BASENAME}"
    
    [[ -f "${REF_DIR}/${GENE_LOC}" ]] || { echo "Missing ${GENE_LOC}"; exit 1; }
    [[ -f "${SUMS_DIR}/${BASENAME}.genes.out" ]] || { echo "Missing ${BASENAME}.genes.out"; exit 1; }

    # Add gene symbols to gene names and save as a new file

    # NR==FNR ensures that we create the map only using the first file
    # For each line of the first file (gene location file), store a mapping from Entrez ID ($1) to gene symbol ($6)
    # The 'next' command skips processing the rest of the code for the first file
    # For the second file (main gene stats file):
    #   - Check if the Entrez ID ($1) exists in the map
    #   - If yes, append the corresponding gene symbol to the end of the line
    #   - If no, append 'NA' to indicate no match was found
    # Finally, print the modified line
    awk 'NR==FNR {gene[$1]=$6; next} {if($1 in gene) $0=$0"\t"gene[$1]; else $0=$0"\tNA"; print}' \
    "${REF_DIR}/${GENE_LOC}" \
    "${SUMS_DIR}/${BASENAME}.genes.out" > "${SUMS_DIR}/${BASENAME}.symbol.out"

    # ---- Generate .gs file -----

    # Build comma-separated GENESET for top 1000 (only positive z-score)
    geneset=$(
        awk 'NR>1 && $10 != "NA" && $8>0 {print $10 ":" $8 "\t" $8}' \
            "${SUMS_DIR}/${BASENAME}.symbol.out" \
        | sort -k2,2nr \
        | cut -f1 \
        | head -n 1000 \
        | paste -sd "," -
    )

    # Skip this trait if no genes were mapped
    [[ -z "$geneset" ]] && continue

    # Append to output file
    echo -e "${BASENAME}\t${geneset}" >> "${OUT_DIR}/all_traits_geneset.gs"
done