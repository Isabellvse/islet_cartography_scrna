#!/bin/bash

# Activate conda environment
# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography_clone

# Load study and necessary paths
study_name="Kang"
cd "/work/islet_cartography_scrna/scripts/download_alignment/${study_name}/"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/scripts/hg38/"
whitelist="/work/islet_cartography_scrna/whitelist/3M-february-2018.txt"
Donors=$(cut -f 1 "$Study" | sort | uniq)

# Load STAR genome into memory
STAR --genomeDir "$Genome" --genomeLoad LoadAndExit
rm Log* Aligned.out.sam SJ.out.tab

# Loop
for z in $Donors; do
	# Setup the files needed to be downloaded for each donor
    awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > Donor
    awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' Donor > Download

     # Download using aria2
    aria2c -i Download -j10 --check-integrity true --save-session failed.downloads
    has_error=`wc -l < failed.downloads`
    
    while [ $has_error -gt 0 ]
    do
    
    echo "still has $has_error errors, rerun aria2 to download ..."
    mv failed.downloads Download
    aria2c -iDownload -j10 --check-integrity true -c --save-session failed.downloads
    has_error=$(wc -l < failed.downloads)
    sleep 10
    done

    for i in *.sra; do
        fasterq-dump "$i" -m 10gb -e 30 -p -S --include-technical
    done

    # Run STAR
    cat *_2.fastq* > R1.fq # barcode + umi
    cat *_3.fastq* > R2.fq # cDNA
    rm *.fastq

    STAR --genomeLoad LoadAndKeep --genomeDir $Genome --readFilesIn R2.fq R1.fq --soloType CB_UMI_Simple --soloFeatures Gene GeneFull GeneFull_Ex50pAS --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --soloCBwhitelist $whitelist --soloBarcodeReadLength 0 --soloUMIlen 12
	
    # Cleanup: Move results to cell-specific folder
    mkdir -p "$Out/$z/"
    mv failed.downloads Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab Download Donor *.fq *.sra
done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab 