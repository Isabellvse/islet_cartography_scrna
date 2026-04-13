#!/usr/bin/env bash
# ----------------------------- INFO -----------------------------------------
# Gene and reference location files were downloaded from: https://cncr.nl/research/magma/ (access: 2025-12-19)
# Installing magma
# Install magma first
# mkdir MAGMA
# cd MAGMA
# wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip
# unzip magma_v1.10.zip
# https://hugeamp.org/help.html?page=955#:~:text=appears%20in%20the%20following%20documents,from%20the%201000%20Genomes%20Project:
# ----------------------------- CONDA -----------------------------------------
eval "$(conda shell.bash hook)"
conda activate "/work/islet_cartography_scrna/scrna_cartography_gwas"

# ----------------------------- USER CONFIG -----------------------------------
# Define folders
WD="/work/islet_cartography_scrna"
TMP_DIR="$WD/tmp"
MAGMA_PROG="$WD/MAGMA"
REF_DIR="$WD/data/gwas/ref_magma"
GWAS_DIR="$WD/data/gwas"
FILE_DIR="$GWAS_DIR/files"
OUT_DIR="$GWAS_DIR/sumstats_annotated"

PROCESS_FILE="$FILE_DIR/test_2.txt"

# Generate folders
mkdir -p "$TMP_DIR"
cd "$TMP_DIR"

# ----------------------------- sumstat files ---------------------------------
files=$(awk -F'\t' 'NR>1 {print $1}' "$PROCESS_FILE" | sort -u)

# ----------------------------- loop ------------------------------------------
for z in $files; do

    echo "Downloading ${z}"
    
    awk -v VAR="$z" -F'\t' '$1 == VAR' "$PROCESS_FILE" > trait

    # --- download file: url, out, checksum ---
    awk -F'\t' '{print $3; print " out=" $2; print " checksum=md5=" $4}' trait > download

    # --- Download ---
    aria2c -i download -j 10 --check-integrity=true --save-session=failed.downloads

    # --- Genome version ---
    genome_version=$(awk -F'\t' '{print $10}' trait | head -1)

    # --- Annotate SNPs ---
    echo "Annotating SNPs"

    # Define N
    read -r N FILE_NAME < <(awk -F'\t' '{print $13, $2}' trait | tr -d '\r\n')
    echo "N = $N"
    echo "FILE_NAME = $FILE_NAME"
    
    python3 <<EOF
import gwaslab as gl

z = "$z"
file_name = "$FILE_NAME"
n = int("$N")

genome_version = "$genome_version"
build = "19" if genome_version == "GRCh37" else "38"

# --- Load file ---
ss = gl.Sumstats(f"{file_name}", n = n, build=build, sep="\t", fmt="ssf", verbose=False)
ss.basic_check(verbose=False)

# Check if SNP ID is there
ss_cols = ss.data.columns
# --- Do liftover from hg38 to hg19 if needed ---
# Only liftover if using GRCh38-style input
if build == "38":
    ss.liftover(from_build="38", to_build="19", remove=True)

# --- Annotate snps ---
ss.assign_rsid(
    ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
    threads=60,
    verbose=False
)

# --- Remove snps with no ID ---
ss.data = ss.data[ss.data['rsID'].notna()]

# --- Rename column rsID -> SNP ---
rename_col = {'rsID': 'SNP'}
ss.data = ss.data.rename(columns=rename_col)

# --- Fill in missing data ---
# If pval or zscore are missing
ss_cols = ss.data.columns

if "P" not in ss.data.columns:
    ss.fill_data(to_fill=["P"])

# --- Save --- 
ss.data.to_csv(f"{z}.pval", sep="\t", index=False)
EOF

    echo "Running MAGMA"
    
    # ---Annotate genes with MAGMA---
    GWAS_PVAL="${z}.pval"
    read -r REF_BIM GENE_LOC < <(awk -F'\t' '{print $8, $9}' trait)

    # - Generate prefix for reference build -
    REF_PREFIX="${REF_BIM%.bim}"
    OUT_PREFIX="${OUT_DIR}/$z"

   echo "REF_BIM = $REF_BIM"
   echo "GENE_LOC = $GENE_LOC"
   echo "REF_PREFIX = $REF_PREFIX"
   echo "GWAS_PVAL = $GWAS_PVAL"
    
    # - Annotation step -
    "${MAGMA_PROG}/magma" --annotate window=10,10 \
        --snp-loc "${REF_DIR}/${REF_BIM}" \
        --gene-loc "${REF_DIR}/${GENE_LOC}" \
        --out "$z"

    # - Gene analysis step -
    "${MAGMA_PROG}/magma" --bfile "${REF_DIR}/${REF_PREFIX}" \
        --pval "${GWAS_PVAL}" use='SNP,P' ncol='N' \
        --gene-annot "$z.genes.annot" \
        --out "${OUT_PREFIX}"

    echo "Add Gene symbol to .out file"
    
    # - Add gene symbol to genes.out file -
    # NR==FNR{gene[$1]=$6;next} → read reference file and build gene ID → symbol map
    # NR==1{print $0"\tSYMBOL";next} → add header for the new column
    # {print $0"\t"(($1 in gene)?gene[$1]:"SYMBOL")} → for each row, append symbol if found, else "NA"

    awk 'NR==FNR{gene[$1]=$6;next} NR==1{print $0"\tSYMBOL";next} {print $0"\t"(($1 in gene)?gene[$1]:"NA")}' \
    "${REF_DIR}/${GENE_LOC}" \
    "${OUT_PREFIX}.genes.out" > "${OUT_PREFIX}.symbol.out"

    echo "Clean up"
    
    # --- Clean up ---
    RESULT_FILES=("${OUT_PREFIX}.genes.out" "${OUT_PREFIX}.symbol.out" "${OUT_PREFIX}.genes.raw" "${OUT_PREFIX}.log.suppl")

    ALL_EXIST=true
    for f in "${RESULT_FILES[@]}"; do
        [ -f "$f" ] || { ALL_EXIST=false; break; }
    done
    
    if $ALL_EXIST; then
        rm -f trait download failed.downloads *.tsv.gz *.genes.annot *.pval
        mv *.log ${OUT_DIR}/
    fi

    echo "Finished ${z}, you can find output files in ${OUT_PREFIX}"

done




