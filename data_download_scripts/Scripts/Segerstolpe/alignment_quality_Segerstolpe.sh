#!/bin/bash

# Activate conda environment
# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

study_name="Segerstolpe"
study="${study_name}.wget"
data="/work/scRNAseq/${study_name}/Preprocessed"

# Output folder
out="/work/islet_cartography_scrna/data_export/${study_name}"
mkdir -p "$out"

# Temporary folder
tmp="${data}/tmp"
mkdir -p "$tmp"

# Create CSV file for features and add prefix gene
awk '{print $1 "_gene_feat" "," $2}' "${data}//Solo.out/Gene/Features.stats" > "$tmp/gene_features.csv"

# Create CSV file for features and add prefix genefull
awk '{print $1 "_genefull_feat" "," $2}' "${data}/Solo.out/GeneFull/Features.stats" > "$tmp/genefull_features.csv"

# Add prefix to summary CSV file gene
awk -F',' '{print $1 "_gene_sum" "," $2}' "${data}/Solo.out/Gene/Summary.csv" > "$tmp/gene_summary.csv"

# Add prefix to summary CSV file genefull
awk -F',' '{print $1 "_genefull_sum" "," $2}' "${data}/Solo.out/GeneFull/Summary.csv" > "$tmp/genefull_summary.csv"

# Create CSV file for barcode.stats
    awk '{print $1 "_barcode" "," $2}' "${data}/${z}/Solo.out/Barcodes.stats" > "$tmp/barcode.csv"

    # Process log file 
    # divide file by | and remove first 4 rows
    awk -F'|' 'NR > 4 {
    # Divide by |
        gsub(/^ +| +$/, "", $1);
        gsub(/^ +| +$/, "", $2);
    # replace white space with _    
        gsub(/ /, "_", $1);
        gsub(/ /, "_", $2);
    # replace , with _
        gsub(/,/, "_", $1);
        gsub(/,/, "_", $2);
    # remove :
        gsub(/:/, "", $1);
        if ($1 != "" && $2 != "") {
            print $1 "," $2
        }
    }' "${data}/${z}/Log.final.out" > "$tmp/log.csv"

    # Create combined file
    echo "sample, $z" > "$tmp/${z}_combined.csv"

    # Set executable permission for the CSV files
    chmod +x "$tmp"/*.csv

    # Combine the files
    cat "$tmp/log.csv" "$tmp/barcode.csv" "$tmp/gene_features.csv" "$tmp/gene_summary.csv" "$tmp/genefull_features.csv" "$tmp/genefull_summary.csv" >> "$tmp/${study_name}_combined.csv"

    # Clean up intermediate files for this donor
    rm "$tmp/log.csv" "$tmp/barcode.csv" "$tmp/gene_features.csv" "$tmp/gene_summary.csv" "$tmp/genefull_features.csv" "$tmp/genefull_summary.csv"

# Transpose
datamash -t, transpose < "$tmp/${study_name}_combined.csv" > "$out/alignment_qc_${study_name}.csv"

# Clean up
rm -r "$tmp"