#!/bin/bash

# Activate conda environment
# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Load study and necessary paths
study_name="Shrestha"
cd "/work/islet_cartography_scrna/scripts/download_alignment/${study_name}/"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Donors=$(cut -f 1 "$Study" | sort | uniq)
Genome="/work/islet_cartography_scrna/scripts/hg38/"
whitelist="/work/islet_cartography_scrna/whitelist/737K-august-2016.txt"

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

	# Unpack the files
	files=$(ls *.sra)
	for i in $files; do fasterq-dump $i -m 10gb -e 20 -p -S --include-technical; done

	# Run STAR
	R1=$(ls *sra_1.fastq)
        R1=$(echo $R1 | tr ' ' ',')
	R2=$(ls *sra_2.fastq)
	R2=$(echo $R2 | tr ' ' ',')
    
	STAR --genomeDir $Genome --readFilesIn $R2 $R1 --soloType CB_UMI_Simple --soloFeatures Gene GeneFull --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --soloCBwhitelist $whitelist --soloBarcodeReadLength 0

	# Cleanup
    mkdir -p "$Out/$z/"
    mv Solo.out Log* "$Out/$z/" 
    rm SJ.out.tab *.fastq Donor Download *.sra Aligned.out.sam failed.download
done
