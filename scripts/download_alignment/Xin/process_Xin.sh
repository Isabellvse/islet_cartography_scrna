#!/bin/bash

# Define variables
study_name="Xin"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/scripts/hg38/"

# Load STAR genome into memory
STAR --genomeDir "$Genome" --genomeLoad LoadAndExit
rm Log* Aligned.out.sam SJ.out.tab

cells=$(cut -f 1 "$Study" | sort | uniq)

for z in $cells; do

# Setup the files needed to be downloaded for each cell

    awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > cell
    awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' cell > Download
    
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

    # Create manifest
    VAR=$(ls *.gz)
    for i in $VAR; do echo $i | sed "s/.fq.gz//g"  >> Cells; done
    cat Cells | sort | uniq -c | awk '{ print $2 }' > ValidCells
    rm Cells
    VAR=$(cat ValidCells)
    for i in $VAR; do echo -e "$i.fq.gz\t-\t$i" > manifest; done
    rm ValidCells

    
    # Run STAR
	STAR --genomeLoad LoadAndKeep --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat


    # Cleanup: Move results to cell-specific folder
    mkdir -p "$Out/$z/"
    mv failed.downloads Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab Download cell *.fq.gz manifest
done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab 