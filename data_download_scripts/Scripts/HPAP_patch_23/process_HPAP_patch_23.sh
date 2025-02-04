#!/bin/bash

# Redirect output and error to a log file
exec > >(tee -i /work/islet_cartography_scrna/data_download_scripts/Scripts/HPAP_patch_23/output.log)
exec 2>&1

cd /work/islet_cartography_scrna/data_download_scripts/Scripts/HPAP_patch_23/

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Define variables
study_name="HPAP_patch_23"
Study="${study_name}.wget"
#Out="/work/scRNAseq/${study_name}/Preprocessed"
#mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"

# check that all files excists
awk '{ print $3}' "$Study" > dry_run
aria2c --dry-run=true -i dry_run --save-session not_available

# extract base name of files not available
grep -o -E 'https?://[^"]+' not_available | awk -F/ '{print $NF}' | awk '{gsub(/-R[12]_fastq-data.fastq.gz/, ""); print}' | awk -F_ '{print $3}' > base_names

# remove files form wget file that do not have all files for each cell (base name)
awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' base_names FS="\t" OFS="\t" "$Study" > "${study_name}_filt.wget"

# define new study
Study="${study_name}_filt.wget"

# clean up
mv not_available base_names $Out
rm dry_run

# Load STAR genome into memory
STAR --genomeDir "$Genome" --genomeLoad LoadAndExit
rm Log* Aligned.out.sam SJ.out.tab

cells=$(cut -f 1 "$Study" | sort | uniq)

for z in $cells; do

# Setup the files needed to be downloaded for each cell
    awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > cell
    awk '{ print $3"\n out="$2 }' cell > Download
    
    # Download using aria2
    aria2c -i Download -j10 --save-session failed.downloads
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
    VAR=$(ls *.gz)
    for i in $VAR; do echo $i | sed "s/_R1.fq.gz//g" - | sed "s/_R2.fq.gz//g" >> Cells; done
    cat Cells | sort | uniq -c | awk '$1 == "2" { print $2 }' - > ValidCells
    rm Cells
    VAR=$(cat ValidCells)
    for i in $VAR; do echo -e $i"_R1.fq.gz\t"$i"_R2.fq.gz\t"$i >> manifest; done
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