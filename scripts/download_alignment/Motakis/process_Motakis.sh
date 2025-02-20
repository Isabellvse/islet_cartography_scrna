#!/bin/bash

# Load study and necessary paths
study_name="Motakis"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"

mkdir -p "$Out"
Donors=$(cut -f 1 "$Study" | sort | uniq)
Genome="/work/islet_cartography_scrna/scripts/hg38/"
whitelist_v2="/work/islet_cartography_scrna/whitelist/737K-august-2016.txt"
whitelist_v3="/work/islet_cartography_scrna/whitelist/3M-february-2018.txt"

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

    # Unpack the files
    for i in *.sra; do
        fasterq-dump "$i" -m 10gb -e 30 -p -S --include-technical
    done

    # Extract donor and chemistry 
    Chem=$(awk -v VAR=$z '$1 == VAR { print $5 }' $Study) 
    Chemistry=$(echo "$Chem" | sort | uniq -c | awk '{ print $2 }')

    # Run STAR
    cat *_2.fastq* > R1.fq
    cat *_3.fastq* > R2.fq
    rm *.fastq

    if [[ "$Chemistry" == "v2" ]]; then
        STAR --genomeLoad LoadAndKeep --genomeDir "$Genome" --readFilesIn R2.fq R1.fq --soloType CB_UMI_Simple --soloFeatures Gene GeneFull \
            --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 \
            --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
            --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --soloCBwhitelist "$whitelist_v2" \
            --soloBarcodeReadLength 0 --soloUMIlen 10
    else
        STAR --genomeLoad LoadAndKeep --genomeDir "$Genome" --readFilesIn R2.fq R1.fq --soloType CB_UMI_Simple --soloFeatures Gene GeneFull \
            --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 \
            --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
            --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --soloCBwhitelist "$whitelist_v3" \
            --soloBarcodeReadLength 0 --soloUMIlen 12
    fi

    # Cleanup: Move results to donor-specific folder
    mkdir -p "$Out/$z/"
    mv failed.downloads Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab *.fastq *.fq Donor Download *.sra
done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab