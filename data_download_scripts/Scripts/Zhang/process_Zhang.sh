#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Load study and necessary paths
study_name="Zhang"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir $Out
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"
Donors=$(cut -f 1 $Study | sort | uniq)

# Loop over donors
for z in $Donors; do
	# Setup the files needed to be downloaded
	awk -v VAR=$z '$1 == VAR { print $0 }' $Study > Donor
	awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' Donor > Download

        # Download using aria2
        aria2c -iDownload -j10 --check-integrity true --save-session failed.downloads
        has_error=`wc -l < failed.downloads`
        while [ $has_error -gt 0 ]
        do
                echo "still has $has_error errors, rerun aria2 to download ..."
                mv failed.downloads Download
                aria2c -iDownload -j10 --check-integrity true -c --save-session failed.downloads
                has_error=`wc -l < out.txt`
                sleep 10
        done

	# Run STAR
	STAR --genomeDir $Genome --readFilesIn *R2.fastq.gz *R1.fastq.gz --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen \ 
    12 --soloUMIstart 13 --soloUMIlen 8 --soloBarcodeReadLength 0 --soloCellFilter None --soloUMIfiltering MultiGeneUMI_CR --runThreadN 20 \
    --outMultimapperOrder Random --outSAMmultNmax 1 --soloCBwhitelist None --readFilesCommand zcat --soloFeatures Gene GeneFull Velocyto \ 
    --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIdedup 1MM_CR 

     # Cleanup: Move results to donor-specific folder
    donor_out_dir="$Out/$z"
    mkdir -p "$donor_out_dir"
    mv Solo.out "$donor_out_dir/"
    mv Log* "$donor_out_dir/"
    rm Aligned.out.sam SJ.out.tab

    # Alignment quality check
    echo "===============================================" >> "$result_file"
    echo "STAR results for $z" >> "$result_file"
    echo "===============================================" >> "$result_file"

    failed_samples=false

    # Process barcode stats 
    barcode_stats="$donor_out_dir/Solo.out/Barcodes.stats"
    if [[ -f $barcode_stats ]]; then
        # Extract barcode stats (nExactMatch and sum of other values)
        valid_barcode_stats=$(awk '
            BEGIN {nExactMatch = 0; sum_others = 0}
            $1 == "nExactMatch" {nExactMatch = $2}
            $1 != "nExactMatch" {sum_others += $2}
            END {
                if (sum_others == 0) {print 0}
                else {print nExactMatch / sum_others}
            }' "$barcode_stats")

        # Check if valid barcode stats pass (must be greater than 1)
        valid_barcode_stats_passed=$(echo "$valid_barcode_stats > 1" | bc)
        if [[ $valid_barcode_stats_passed -eq 1 ]]; then
            echo "Barcodes stats passed with nExactMatch / sum_others: $valid_barcode_stats" >> "$result_file"
        else
            echo "Invalid barcode stats: nExactMatch / sum_others: $valid_barcode_stats" >> "$result_file"
            failed_samples=true
        fi
    else
        echo "$z: Barcodes.stats not found!" >> "$result_file"
        failed_samples=true
    fi

    # Loop through each folder for STAR results
    for folder in "${folders[@]}"; do
        folder_path="$donor_out_dir/Solo.out/$folder"
        summary_file="$folder_path/Summary.csv"

        # Check if the summary.csv exists and process barcode stats
        if [[ -f $summary_file ]]; then
            # Extract barcode stats
            valid_barcodes=$(awk -F',' '/Reads With Valid Barcodes/ {print $2}' "$summary_file")
            mapped_to_genome=$(awk -F',' '$1 == "Reads Mapped to Genome: Unique" {print $2}' "$summary_file")

            # Check if values pass conditions (using bc for floating point comparison)
            valid_barcodes_passed=$(echo "$valid_barcodes >= 0.8" | bc)
            mapped_to_genome_passed=$(echo "$mapped_to_genome >= 0.8" | bc)

            if [[ $valid_barcodes_passed -eq 1 && $mapped_to_genome_passed -eq 1 ]]; then
                echo "$folder passed with Reads With Valid Barcodes: $valid_barcodes and Reads Mapped to Genome: $mapped_to_genome" >> "$result_file"
            else
                echo "$folder failed. Reads With Valid Barcodes: $valid_barcodes, Reads Mapped to Genome: $mapped_to_genome" >> "$result_file"
                failed_samples=true
            fi
        else
            echo "$folder: summary.csv not found!" >> "$result_file"
            failed_samples=true
        fi
    done

    # Cleanup: Only remove .fastq files if no folder failed
    if [ "$failed_samples" = false ]; then
        echo "No failures detected, proceeding with cleanup..."
        rm *.fastq
    else
        echo "Samples failed, skipping removal of .fastq files"
    fi

    # Final cleanup for this donor
    rm Donor Download *.sra
done