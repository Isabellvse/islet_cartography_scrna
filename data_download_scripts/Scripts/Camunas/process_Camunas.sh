#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Define variables
study_name="Camunas"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"

Donors=$(cut -f 1 "$Study" | sort | uniq)

for z in $Donors; do
    # Setup the files needed to be downloaded for each donor
    awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > Donor
    awk '{  print $4"\n out="$3"\n checksum=md5="$5 }' Donor > Download

     # Download using aria2
    aria2c -i Download -j10 --check-integrity=true --save-session failed.downloads
    has_error=$(wc -l < failed.downloads)
    
   # while [ $has_error -gt 0 ]; do
        echo "Still has $has_error errors, rerun aria2 to download ..."
        mv failed.downloads Download
        aria2c -i Download -j10 --check-integrity=true -c --save-session failed.downloads
        has_error=$(wc -l < failed.downloads)
        sleep 10
   # done

    # Create a manifest file
    VAR=$(ls *.fq.gz)
    for i in $VAR; do cell=$(echo $i | sed 's/.fq.gz//g'); echo -e $i"\t-\t"$cell >> manifest; done

    # run STAR 
    STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

# Cleanup
mv Solo.out $Out
mv Log* $Out
mv failed.downloads $Out

rm Aligned.out.sam
rm SJ.out.tab
rm *.fq.gz*
rm Donor
rm Download
rm manifest

done 