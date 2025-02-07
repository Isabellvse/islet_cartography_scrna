#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Define variables
study_name="Baron"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"
bc1="/work/islet_cartography_scrna/whitelist/indrop_bc1.txt"
bc2="/work/islet_cartography_scrna/whitelist/indrop_bc2.txt"

Donors=$(cut -f 1 "$Study" | sort | uniq)


# Load STAR genome into memory
STAR --genomeDir "$Genome" --genomeLoad LoadAndExit
rm Log* Aligned.out.sam SJ.out.tab


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
		
	# Combine multiple runs
	cat *_R1.fq.gz > R1.fq.gz
	cat *_R2.fq.gz > R2.fq.gz
	
	# Run STAR
	STAR --genomeLoad LoadAndKeep --genomeDir $Genome --readFilesIn R2.fq.gz R1.fq.gz --soloType CB_UMI_Complex --soloCBwhitelist $bc1 $bc2 --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT --soloAdapterMismatchesNmax 2 --soloCBposition 0_0_2_-1 3_1_3_8 --soloUMIposition 3_9_3_14 --soloCBmatchWLtype 1MM --soloBarcodeReadLength 0 --soloCellFilter None --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloFeatures Gene GeneFull Velocyto --outFilterScoreMin 30 --soloMultiMappers EM --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

	# Cleanup
	mkdir -p "$Out/$z/"
    mv failed.downloads Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab Download Donor *.fq.gz 
done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab


