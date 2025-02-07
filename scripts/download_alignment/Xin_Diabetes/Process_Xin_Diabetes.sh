#!/bin/bash
# Load information about the number of samples in the study
Out="/data/home/jgsm/VDS2/Mandrup/JM/Pancreas/Data/scRNAseq/Xin_Diabetes/"
Study="Xin_Diabetes.wget"
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

	# Create Fastqs from 10X BAM files
	samtools view -F 256 *.bam | awk '{ print "@"$1" 1:N:0:0\n"$10"\n+\n"$11 }' > R2.fq
	samtools view -F 256 -x GX -x GN -x RE -x CB -x BC -x QT -x RG -x UB -x MM *.bam | awk '{ print "@"$1" 2:N:0:0\n"$(NF-3)$(NF-1)"\n+\n"$(NF-2)$(NF) }' | sed 's/CR:Z://g' | sed 's/UR:Z://g' | sed 's/CY:Z://g' | sed 's/UY:Z://g' > R1.fq

	# Run STAR
	STAR --genomeDir $Genome --readFilesIn R2.fq R1.fq --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10 --soloBarcodeReadLength 1 --soloCellFilter None --soloUMIfiltering MultiGeneUMI_CR --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1  --soloCBwhitelist Xin.whitelist --soloUMIdedup 1MM_CR --soloFeatures Gene GeneFull Velocyto --outFilterScoreMin 30 --soloMultiMappers EM --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4

	# Cleanup
	mkdir $Out/$z/
	mv Solo.out $Out/$z/
	rm Aligned.out.sam
	mv Log* $Out/$z/
	rm SJ.out.tab
	rm *.fq*
	rm Donor
	rm *.bam*
done
