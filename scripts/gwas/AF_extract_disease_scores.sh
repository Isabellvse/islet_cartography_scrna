#!/usr/bin/env bash
# ----------------------------- CONDA -----------------------------------------
eval "$(conda shell.bash hook)"
conda activate /work/islet_cartography_scrna/scrna_cartography_gwas

# ----------------------------- USER CONFIG -----------------------------------
WD="/work/islet_cartography_scrna"
GWAS_DIR="$WD/data/gwas"
OUT_DIR="$GWAS_DIR/disease_scores"
TMP_DIR="$WD/tmp"

mkdir -p "$TMP_DIR"
cd "$TMP_DIR"

# ----------------------------- CREATE TRAITS FILE -----------------------------------
DS_DIR="$OUT_DIR"  # directory with <trait>_pops.full_score.gz files
TRAITS_FILE="$TMP_DIR/traits.txt"

> "$TRAITS_FILE"
shopt -s nullglob
for f in "$DS_DIR"/*_pops.full_score.gz; do
    trait=$(basename "$f" "_pops.full_score.gz")
    echo "$trait" >> "$TRAITS_FILE"
done
shopt -u nullglob

# ----------------------------- PREPARE TMP DIRECTORIES -----------------------------------
mkdir -p "$TMP_DIR/aligned_scores"

export DS_DIR TMP_DIR

# ----------------------------- EXTRACT AND SORT SCORES BY BARCODE -----------------------------------
xargs -a "$TRAITS_FILE" -P 8 -I{} sh -c '
zcat "$DS_DIR/{}_pops.full_score.gz" | awk -F"\t" "
NR==1 {
  for(i=1;i<=NF;i++){ if(\$i==\"norm_score\") ncol=i; if(\$i==\"raw_score\") rcol=i }
}
NR>1 {
  print \$1\"\t\"\$ncol\"\t\"\$rcol
}
" | sort -k1,1 > "$TMP_DIR/aligned_scores/{}.txt"
'

# ----------------------------- MERGE FILES BY BARCODE -----------------------------------
# Take the first trait as the base
read -r first_trait < "$TRAITS_FILE"
cp "$TMP_DIR/aligned_scores/$first_trait.txt" "$TMP_DIR/merged.txt"

# Loop over remaining traits and join by barcode
# sort files based on column 1
tail -n +2 "$TRAITS_FILE" | while read -r trait; do
    join -t $'\t' "$TMP_DIR/merged.txt" "$TMP_DIR/aligned_scores/$trait.txt" > "$TMP_DIR/merged_tmp.txt"
    mv "$TMP_DIR/merged_tmp.txt" "$TMP_DIR/merged.txt"
done

# ----------------------------- CREATE HEADER -----------------------------------
# Save header
header="barcode"
while read -r trait; do
    header+=$'\t'"${trait}_norm"$'\t'"${trait}_raw"
done < "$TRAITS_FILE"

# ----------------------------- SAVE FINAL MERGED FILE -----------------------------------
echo -e "$header" | cat - "$TMP_DIR/merged.txt" > "$OUT_DIR/merged_scores.txt"