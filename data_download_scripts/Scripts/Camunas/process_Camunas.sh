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

# Process data
awk '{  print $3"\n out="$2"\n checksum=md5="$4 }' $Study > Download

# Download using aria2
aria2c -i Download -j10 --check-integrity true --save-session failed.downloads
has_error=`wc -l < failed.downloads`
while [ $has_error -gt 0 ]
do
echo "still has $has_error errors, rerun aria2 to download ..."
mv failed.downloads Download
aria2c -iDownload -j10 --check-integrity true -c --save-session failed.downloads
has_error=`wc -l < out.txt`
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
STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 50 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

# Cleanup
mv Solo.out $Out
rm Aligned.out.sam
mv Log* $Out
rm SJ.out.tab
rm *.fq.gz*
rm Donor
rm Download
rm manifest
mv failed.downloads $Out