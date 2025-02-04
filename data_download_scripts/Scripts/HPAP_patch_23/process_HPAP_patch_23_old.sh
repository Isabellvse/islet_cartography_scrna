#!/bin/bash

# Internalize shell
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate /work/islet_cartography_scrna/scrna_cartography

# Define variables
study_name="HPAP_patch_23"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"

# check that all files excists
awk '{ print $3}' "$Study" > dry_run
aria2c --dry-run=true -i dry_run --save-session not_available

# extract base name of files not available
grep -o -E 'https?://[^"]+' not_available | awk -F/ '{print $NF}' | awk '{gsub(/-R[12]_fastq-data.fastq.gz/, ""); print}' > base_names

# remove files form wget file that do not have all files for each cell (base name)
awk 'NR==FNR {exclude[$1]; next} {base=substr($3, 1, index($3, "-R")-1)} !(base in exclude)' base_names FS="\t" OFS="\t" "$Study" > "${study_name}_filt.wget"

# define new study
Study="${study_name}_filt.wget"

# clean up
mv not_available $Out
mv base_names $Out
rm dry_run

# Extract donors
Donors=$(cut -f 1 "$Study" | sort | uniq)

for z in $Donors; do
  echo "Processing donor: $z"
  
  # Process data
  awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > Donor
  awk '{ print $4 "\n out=" $2".fq.gz" }' Donor > Download

    # Download using aria2 with retry logic - keep trying to download untill all files are downloaded successfully
  aria2c -i Download -j10 --save-session failed.downloads --on-download-error="echo 'Download error for $z' >> download_errors.log"
  has_error=$(wc -l < failed.downloads)
  while [ $has_error -gt 0 ]; do
    echo "still has $has_error errors, rerun aria2 to download ..."
    mv failed.downloads Download
    aria2c -i Download -j10 -c --save-session failed.downloads --on-download-error="echo 'Download error for $z' >> download_errors.log"
    has_error=$(wc -l < failed.downloads)
    sleep 10
  done
  
  # Create the manifest file
  VAR=$(ls *.gz)
  for i in $VAR; do 
    echo $i | sed "s/_R1.fq.gz//g" - | sed "s/_R2.fq.gz//g" >> Cells
  done
  cat Cells | sort | uniq -c | awk '$1 == "2" { print $2 }' - > ValidCells
  rm Cells
  VAR=$(cat ValidCells)
  for i in $VAR; do 
    echo -e $i"_R1.fq.gz\t"$i"_R2.fq.gz\t"$i >> manifest
  done
  rm ValidCells
  
  echo "Manifest file created for donor: $z"
  
  # Print first part of manifest
  head manifest

  # Run STAR
  STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat
  
  echo "STAR alignment completed for donor: $z"
  
  # Move the results to correct folder
  mkdir $Out/$z
  mv Solo.out/ $Out/$z
  mv Log.* $Out/$z
  mv Err.file $Out/$z
  mv failed.downloads $Out/$z

  echo "Results moved for donor: $z"

  # Cleanup the processing folder
  rm *.gz
  rm Aligned.out.sam
  rm SJ.out.tab
  rm manifest
  rm Donor
  rm Download

  echo "Cleanup completed for donor: $z"
done