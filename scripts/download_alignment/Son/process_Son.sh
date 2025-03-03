#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography
exec > >(tee -i /work/islet_cartography_scrna/scripts/download_alignment/Son/output.log)
exec 2>&1
cd /work/islet_cartography_scrna/scripts/download_alignment/Son/

# Define variables
study_name="Son"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/scripts/hg38/"

# Dry run
#awk '{ print $3}' "$Study" > dry_run
#aria2c --dry-run=true -i dry_run --check-integrity true --save-session not_available

# Get urls that were not available
#grep -o 'ftp://[^ ]*' not_available > urls_to_remove.txt

# Copy Son.wget to Son_filt.wget
#cp Son.wget Son_filt.wget

# Remove lines from Son_filt.wget
#while read -r url; do
#    grep -v "$url" Son_filt.wget > temp && mv temp Son_filt.wget
#done < urls_to_remove.txt

# Clean up
#rm dry_run
#mv not_available urls_to_remove.txt "$Out/"

# Define new wget file
Study="${study_name}_filt.wget"

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

    # Create the manifest file
    VAR=$(ls *.fq.gz)
    for i in $VAR; do echo $i | sed "s/.fq.gz//g" > Cells; done
    VAR=$(cat Cells)
    for i in $VAR; do echo -e "$i".fq.gz"\t-\t$i" > manifest; done
    rm Cells
    
    # Run STAR
    STAR --genomeLoad LoadAndKeep --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Forward --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

    # Cleanup: Move results to cell-specific folder
    mkdir -p "$Out/$z/"
    mv Solo.out Log* "$Out/$z/"
    rm Aligned.out.sam SJ.out.tab Download cell *.fq.gz manifest

done

# Unload the genome from shared memory
STAR --genomeDir "$Genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab