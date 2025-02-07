#!/bin/bash

# Load study and necessary paths
study_name="Wang_Sander"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"

mkdir -p "$Out"
Donors=$(cut -f 1 "$Study" | sort | uniq)
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"
whitelist="/work/islet_cartography_scrna/whitelist/737K-arc-v1.txt"

# Load STAR genome into memory
STAR --genomeDir "$Genome" --genomeLoad LoadAndExit
rm Log* Aligned.out.sam SJ.out.tab

# Loop over donors
for z in $Donors; do
    echo "Processing donor: $z"

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

    # Run STAR
    STAR --genomeLoad LoadAndKeep --genomeDir "$Genome" --readFilesIn *R2.fq.gz *R1.fq.gz --soloType CB_UMI_Simple --soloFeatures Gene GeneFull --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --soloCBwhitelist "$whitelist" --soloBarcodeReadLength 0 --soloUMIlen 12 --readFilesCommand zcat

    # Cleanup: Move results to donor-specific folder
    mkdir -p "$Out/$z/"
    mv failed.downloads Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab *.fq.gz Donor Download 
done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab