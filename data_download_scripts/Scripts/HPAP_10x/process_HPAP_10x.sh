#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Load study and necessary paths
study_name="HPAP_10x"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"

Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"
whitelist_v2="/work/islet_cartography_scrna/whitelist/737K-august-2016.txt"
whitelist_v3="/work/islet_cartography_scrna/whitelist/3M-february-2018.txt"
Donors=$(cut -f 1 "$Study" | sort | uniq)

# Loop
for z in $Donors; do
	# Setup the files needed to be downloaded for each donor
    awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > Donor
    awk '{ print $3 "\n out=" $2 }' Donor > Download

	# Download the data
    aria2c -i Download -j10

    # Extract donor and chemistry information from the samples files (subset from official HPAP metadata)
	Chemistry=$(awk -v VAR=$z '$1 == VAR { print $4 }' $Study)

	# Cat the individual data files and remove the individual files
	cat *R1* > R1.fq.gz
	cat *R2* > R2.fq.gz
	rm *.fastq.gz
	
	# Align the data using STARsolo
	if [[ "$Chemistry" == "v2" ]];
	then
        	STAR --genomeDir $Genome --readFilesIn R2.fq.gz R1.fq.gz --soloType CB_UMI_Simple --soloFeatures Gene GeneFull --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat --soloCBwhitelist $whitelist_v2 --soloBarcodeReadLength 0
	else
        	STAR --genomeDir $Genome --readFilesIn R2.fq.gz R1.fq.gz --soloType CB_UMI_Simple --soloFeatures Gene GeneFull --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat --soloCBwhitelist $whitelist_v3 --soloBarcodeReadLength 0 --soloUMIlen 12
	fi

    # Move the results
    mkdir $Out/$z
    mv Solo.out/ $Out/$z
    mv Log.* $Out/$z

	# Cleanup the processing folder
	rm *.fq.gz Aligned.out.sam SJ.out.tab
done
