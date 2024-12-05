#!/bin/bash

# Activate conda environment
# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

study_name="Camunas"
study="${study_name}.wget"
data="/work/scRNAseq/${study_name}/Preprocessed"

# Output folder
out="/work/islet_cartography_scrna/data_export/${study_name}"
mkdir -p "$out"

# Temporary folder
tmp="${data}/tmp"
mkdir -p "$tmp"

donors=$(cut -f 1 "$study" | sort | uniq)

for z in $donors; do
    # Create CSV file for features and add prefix gene
    awk '{print $1 "_gene_feat" "," $2}' "${data}/${z}/Solo.out/Gene/Features.stats" > "$tmp/gene_features.csv"

    # Create CSV file for features and add prefix genefull
    awk '{print $1 "_genefull_feat" "," $2}' "${data}/${z}/Solo.out/GeneFull/Features.stats" > "$tmp/genefull_features.csv"

    # Add prefix to summary CSV file gene
    awk -F',' '{print $1 "_gene_sum" "," $2}' "${data}/${z}/Solo.out/Gene/Summary.csv" > "$tmp/gene_summary.csv"

    # Add prefix to summary CSV file genefull
    awk -F',' '{print $1 "_genefull_sum" "," $2}' "${data}/${z}/Solo.out/GeneFull/Summary.csv" > "$tmp/genefull_summary.csv"

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
    cat "$tmp/log.csv" "$tmp/barcode.csv" "$tmp/gene_features.csv" "$tmp/gene_summary.csv" "$tmp/genefull_features.csv" "$tmp/genefull_summary.csv" >> "$tmp/${z}_combined.csv"

    # Clean up intermediate files for this donor
    rm "$tmp/log.csv" "$tmp/barcode.csv" "$tmp/gene_features.csv" "$tmp/gene_summary.csv" "$tmp/genefull_features.csv" "$tmp/genefull_summary.csv"
done

# Combine all files into one
files=($tmp/*combined.csv)

# Extract the header from the first file
head -n 1 "${files[0]}" > "$tmp/header0.csv"

# Extract headers from all files and combine them
for ((i = 1; i < ${#files[@]}; i++)); do
    head -n 1 "${files[i]}" | cut -d, -f2- > "$tmp/header${i}.csv"
done

# Combine headers into a single line
paste -d, "$tmp/header0.csv" $(for ((i = 1; i < ${#files[@]}; i++)); do echo "$tmp/header${i}.csv"; done) > "$tmp/combined_header.csv"

# Remove headers from the original files and sort them
for file in "${files[@]}"; do
    tail -n +2 "$file" | sort > "${file}_sorted"
done

# Initialize the combined body file with the first sorted file
cp "${files[0]}_sorted" "$tmp/combined_body.csv"

# Loop through the rest of the sorted files and join them
for ((i = 1; i < ${#files[@]}; i++)); do
    join -t, -1 1 -2 1 "$tmp/combined_body.csv" "${files[i]}_sorted" > "$tmp/temp_combined.csv"
    mv "$tmp/temp_combined.csv" "$tmp/combined_body.csv"
done

# Combine the header and body
cat "$tmp/combined_header.csv" "$tmp/combined_body.csv" > "$tmp/combined_t.csv"

# Transpose
datamash -t, transpose < "$tmp/combined_t.csv" > "$out/alignment_qc_${study_name}.csv"

# Clean up
rm -r "$tmp"