#!/bin/bash

# Activate conda environment
# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Load study and necessary paths
study_name="Muraro"
cd "/work/islet_cartography_scrna/scripts/download_alignment/${study_name}/"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Donors=$(cut -f 1 "$Study" | sort | uniq)
Genome="/work/islet_cartography_scrna/scripts/hg38/"
whitelist="/work/islet_cartography_scrna/whitelist/Muraro.whitelist"

# Load STAR genome into memory
STAR --genomeDir "$Genome" --genomeLoad LoadAndExit
rm Log*

# Loop over donors
for z in $Donors; do
    # Setup the files needed to be downloaded for each donor
    awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > Donor
    awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' Donor > Download

    # Download using aria2
    aria2c -i Download -j10 --check-integrity=true --save-session failed.downloads
    has_error=$(wc -l < failed.downloads)
    
    while [ $has_error -gt 0 ]; do
        echo "Still has $has_error errors, rerun aria2 to download ..."
        mv failed.downloads Download
        aria2c -i Download -j10 --check-integrity=true -c --save-session failed.downloads
        has_error=$(wc -l < failed.downloads)
        sleep 10
    done
    
    for i in *.sra; do
        fasterq-dump "$i" -m 10gb -e 30 -p -S --include-technical
    done

    # Run STAR
    	STAR --genomeLoad LoadAndKeep --genomeDir $Genome --readFilesIn *_2.fastq *_1.fastq --soloType CB_UMI_Simple --soloCBwhitelist $whitelist --soloCBstart 1 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen 4 --soloBarcodeReadLength 0 --soloCellFilter None --soloUMIfiltering - --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --soloUMIdedup Exact --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCBmatchWLtype 1MM

    # Cleanup: Move results to donor-specific folder
    mkdir -p "$Out/$z/"
    mv failed.downloads Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab Download Donor *.fastq *.sra
done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab