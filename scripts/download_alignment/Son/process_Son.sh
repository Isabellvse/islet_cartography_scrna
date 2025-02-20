#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography
exec > >(tee -i /work/islet_cartography_scrna/scripts/download_alignment/Son/output_sra.log)
exec 2>&1
cd /work/islet_cartography_scrna/scripts/download_alignment/Son/

# Define variables
study_name="Son"
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
    attempt=0
    max_attempts=5
    aria2c -i Download -j10 --check-integrity true --save-session failed.downloads
    has_error=$(wc -l < failed.downloads)

    # Retry logic for failed downloads
    while [ $has_error -gt 0 ] && [ $attempt -lt $max_attempts ]; do
        echo "Attempt $((attempt + 1)) of $max_attempts: still has $has_error errors, rerunning aria2 to download ..."
        mv failed.downloads Download
        aria2c -i Download -j10 --check-integrity true -c --save-session failed.downloads
        has_error=$(wc -l < failed.downloads)
        attempt=$((attempt + 1))
        sleep 10
    done

    # If there are still errors after 5 attempts, skip this cell and move failed.downloads
    if [ $has_error -gt 0 ]; then
        echo "Download failed after $max_attempts attempts, skipping file processing and moving failed.downloads."
        mkdir -p "$Out/$z/"
        mv failed.downloads "$Out/$z/"
        rm Download cell
        continue
    fi

    # Process the downloaded files
    for i in *.sra; do
        fasterq-dump "$i" -m 10gb -e 60 -p -s --include-technical
    done

    rm *sra
    
    # Create the manifest file
    VAR=$(ls *.gz)
    for i in $VAR; do echo $i | sed "s/.fastq//g" > Cells; done
    VAR=$(cat Cells)
    for i in $VAR; do echo -e "$i.fastq\t-\t$i" > manifest; done
    rm Cells
    
    # Run STAR
    STAR --genomeLoad LoadAndKeep --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Forward --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1

    # Cleanup: Move results to cell-specific folder
    mkdir -p "$Out/$z/"
    mv Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab Download cell *fastq manifest

done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab
