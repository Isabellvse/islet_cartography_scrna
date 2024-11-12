#!/bin/bash
# Load information about the number of samples in the study
Out="/data/home/jgsm/VDS2/Mandrup/JM/Pancreas/Data/scRNAseq/Baron/"
Study="Baron.wget"
Donors=$(cut -f 1 $Study | sort | uniq)
Genome="/data/home/jgsm/Pancreas/hg38/"
mkdir $Out

# Loop over donors
for z in $Donors; do
	# Setup the files needed to be downloaded
	awk -v VAR=$z '$1 == VAR { print $0 }' $Study > Donor
	awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' Donor > Download

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

	# Run STAR
	STAR --genomeDir $Genome --readFilesIn *R2.fq.gz *R1.fq.gz --soloType CB_UMI_Complex --soloCBwhitelist bc1.txt bc2.txt --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT --soloAdapterMismatchesNmax 2 --soloCBposition 0_0_2_-1 3_1_3_8 --soloUMIposition 3_9_3_14 --soloCBmatchWLtype 1MM --soloBarcodeReadLength 0 --soloCellFilter None --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloFeatures Gene GeneFull Velocyto --outFilterScoreMin 30 --soloMultiMappers EM --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

	# Cleanup
	mkdir $Out/$z/
	mv Solo.out $Out/$z/
	rm Aligned.out.sam
	mv Log* $Out/$z/
	rm SJ.out.tab
	rm *.fq.gz*
	rm Download
	rm Donor
done
