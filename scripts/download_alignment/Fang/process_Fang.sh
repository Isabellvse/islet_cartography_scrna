#!/bin/bash

# Activate conda environment
# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Load information about the number of samples in the study
study_name="Fang"
cd "/work/islet_cartography_scrna/scripts/download_alignment/${study_name}/"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Donors=$(cut -f 1 $Study | sort | uniq)
Genome="/work/islet_cartography_scrna/scripts/hg38/"

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
	STAR --genomeDir $Genome --readFilesIn *R2.fq.gz *R1.fq.gz --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 --soloBarcodeReadLength 0 --soloCellFilter None --soloUMIfiltering MultiGeneUMI_CR --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --soloCBwhitelist None --readFilesCommand zcat --soloFeatures Gene GeneFull Velocyto --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIdedup 1MM_CR 

	# Cleanup
    mkdir -p $Out/$z/
	mv Solo.out Log* failed.downloads $Out/$z/
	rm Aligned.out.sam SJ.out.tab *.fq.gz* Donor

done
