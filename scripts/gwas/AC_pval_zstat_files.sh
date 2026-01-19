#!/usr/bin/env bash
# ----------------------------- INFO -----------------------------------------


# ----------------------------- CONDA -----------------------------------------
eval "$(conda shell.bash hook)"
conda activate "/work/islet_cartography_scrna/scrna_cartography_gwas"

# ----------------------------- USER CONFIG -----------------------------------
# Define folders
WD="/work/islet_cartography_scrna"
TMP_DIR="$WD/tmp"
GWAS_DIR="$WD/data/gwas"
FILE_DIR="$GWAS_DIR/files"
SUM_FILES="$GWAS_DIR/sumstats_annotated"
PZ_DIR="$GWAS_DIR/pval_zscore"
GS_DIR="$GWAS_DIR/gs_files"

mkdir -p "$TMP_DIR"
cd "$TMP_DIR"
# ----------------------------- LOOPS -----------------------------------------
# Get symbol files
FILES="$SUM_FILES/*symbol.out"
# ----------------------------- PVAL FILES ------------------------------------
for z in $FILES; do
    file_name=$(basename "$z" .symbol.out)

    awk -v col="$file_name" '
        NR==1 {
            for (i=1; i<=NF; i++) if ($i=="P") p=i
            printf "GENE\t%s\n", col
            next
        }
        {
            printf "%s\t%s\n", $NF, $p
        }
    ' "$z" > "${TMP_DIR}/${file_name}_pval.tsv"
done

python3 - <<EOF
import pandas as pd
import glob, os

out_dir="$PZ_DIR"
dfs = []
for f in glob.glob("*pval.tsv"):
    df = pd.read_csv(f, sep="\t")
    dfs.append(df.set_index("GENE"))

out = pd.concat(dfs, axis=1, join="outer")
out.to_csv(os.path.join(out_dir, "pval.tsv"), sep="\t", na_rep="NA")
EOF

# ----------------------------- ZSTAT ------------------------------------
for z in $FILES; do
    file_name=$(basename "$z" .symbol.out)

    awk -v col="$file_name" '
        NR==1 {
            for (i=1; i<=NF; i++) if ($i=="ZSTAT") p=i
            printf "GENE\t%s\n", col
            next
        }
        {
            printf "%s\t%s\n", $NF, $p
        }
    ' "$z" > "${TMP_DIR}/${file_name}_zstat.tsv"
done

python3 - <<EOF
import pandas as pd
import glob, os

out_dir="$PZ_DIR"
dfs = []
for f in glob.glob("*zstat.tsv"):
    df = pd.read_csv(f, sep="\t")
    dfs.append(df.set_index("GENE"))

out = pd.concat(dfs, axis=1, join="outer")
out.to_csv(os.path.join(out_dir, "zstat.tsv"), sep="\t", na_rep="NA")
EOF

# ----------------------------- CLEANUP ------------------------------------
RESULT_FILES=("${PZ_DIR}/pval.tsv" "${PZ_DIR}/zstat.tsv")

ALL_EXIST=true
for f in "${RESULT_FILES[@]}"; do
    [ -f "$f" ] || { ALL_EXIST=false; break; }
done

if $ALL_EXIST; then
    rm -f *.tsv
    echo "Finished generating pval and zstat files, you can find output files in ${PZ_DIR}"
fi

