#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

study_name="Son"
cd "/work/islet_cartography_scrna/scripts/download_alignment/${study_name}/"
study="${study_name}_filt.wget"
data="/work/scRNAseq/${study_name}/Preprocessed"

# Output folder
out="/work/islet_cartography_scrna/data/star_quality"
mkdir -p "$out"

# Temporary folder
tmp="${data}/tmp"
mkdir -p "$tmp"

donors=$(cut -f 1 "$study" | sort | uniq)

process_donor() {
    z=$1

    # Create CSV file for features and add prefix gene
    awk '{print $1 "_gene_feat" "," $2}' "${data}/${z}/Solo.out/Gene/Features.stats" > "$tmp/${z}_gene_features.csv"

    # Create CSV file for features and add prefix genefull
    awk '{print $1 "_genefull_feat" "," $2}' "${data}/${z}/Solo.out/GeneFull/Features.stats" > "$tmp/${z}_genefull_features.csv"

    # Add prefix to summary CSV file gene
    awk -F',' '{print $1 "_gene_sum" "," $2}' "${data}/${z}/Solo.out/Gene/Summary.csv" > "$tmp/${z}_gene_summary.csv"

    # Add prefix to summary CSV file genefull
    awk -F',' '{print $1 "_genefull_sum" "," $2}' "${data}/${z}/Solo.out/GeneFull/Summary.csv" > "$tmp/${z}_genefull_summary.csv"

    # Create CSV file for barcode.stats
    awk '{print $1 "_barcode" "," $2}' "${data}/${z}/Solo.out/Barcodes.stats" > "$tmp/${z}_barcode.csv"

    # Process log file 
    awk -F'|' 'NR > 4 {
        gsub(/^ +| +$/, "", $1);
        gsub(/^ +| +$/, "", $2);
        gsub(/ /, "_", $1);
        gsub(/ /, "_", $2);
        gsub(/,/, "_", $1);
        gsub(/,/, "_", $2);
        gsub(/:/, "", $1);
        if ($1 != "" && $2 != "") {
            print $1 "," $2
        }
    }' "${data}/${z}/Log.final.out" > "$tmp/${z}_log.csv"

    # Create combined file
    echo "sample, $z" > "$tmp/${z}_combined.csv"

    # Combine the files
    cat "$tmp/${z}_log.csv" "$tmp/${z}_barcode.csv" "$tmp/${z}_gene_features.csv" "$tmp/${z}_gene_summary.csv" "$tmp/${z}_genefull_features.csv" "$tmp/${z}_genefull_summary.csv" >> "$tmp/${z}_combined.csv"

    # Clean up intermediate files for this donor
    rm "$tmp/${z}_log.csv" "$tmp/${z}_barcode.csv" "$tmp/${z}_gene_features.csv" "$tmp/${z}_gene_summary.csv" "$tmp/${z}_genefull_features.csv" "$tmp/${z}_genefull_summary.csv"
}

export -f process_donor  # Export function so parallel can use it
export tmp data  # Export necessary variables

# Parallelize the donor processing loop using parallel
echo "$donors" | parallel -j 30 process_donor  # Adjust -j for the number of CPU cores available

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