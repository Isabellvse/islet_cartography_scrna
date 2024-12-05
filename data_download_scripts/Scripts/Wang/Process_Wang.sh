#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Load study and necessary paths
study_name="Wang"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"

# Setup the files needed to be downloaded
awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' $Study > Download

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

# Check for errors
Files=$(ls *.gz)
touch Err.file
for i in $Files; do gzip -t $i 2>&1 | awk 'NR > 1 { print $0 }' - | awk -F ":" '{ print $2 }' - | sed 's/ //g' >> Err.file; done
Del=$(cat Err.file)
for i in $Del; do rm $i; done

# Create a manifest file
VAR=$(ls *.fq.gz)
for i in $VAR; do cell=$(echo $i | sed 's/.fq.gz//g'); echo -e $i"\t-\t"$cell >> manifest; done

# Run STAR
STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

# Cleanup
mv Solo.out $Out
rm Aligned.out.sam
mv Log* $Out
rm SJ.out.tab
rm *.gz*
rm Donor
rm Download
rm manifest
mv failed.downloads $Out
mv Err.file $Out
