#!/bin/bash
# Load information about the number of samples in the study
mkdir /data/home/jgsm/VDS2/Mandrup/JM/Pancreas/Data/scRNAseq/Muraro/
Study="Muraro.txt"
Donors=$(cut -f 1 $Study | sort | uniq)

# Loop over donors
for z in $Donors; do
	# Setup the files needed to be downloaded
	awk -v VAR=$z '$1 == VAR { print $0 }' $Study > Donor		

	# Download using ascp (few big files)
	Files=$(wc -l Donor | awk '{ print $1 }')
        for m in $(seq 1 $Files); do
		Source=$(awk -v VAR=$m 'NR == VAR { print $3 }' Donor)
		Target=$(awk -v VAR=$m 'NR == VAR { print $2 }' Donor)
		ChkSum=$(awk -v VAR=$m 'NR == VAR { print $4 }' Donor)
		NewSum=""
		while [[ "$ChkSum" != "$NewSum" ]]
		do
			/data/home/jgsm/.aspera/connect/bin/ascp -QT -l 300m -m 50m -P33001 -i /data/home/jgsm/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@$Source $Target
			NewSum=$(md5sum $Target | awk '{ print $1 }' - )		
		done
	done

	# Run STAR
	STAR --genomeDir /data/home/jgsm/Pancreas/hg38/ --readFilesIn *R2.fq.gz *R1.fq.gz --soloType CB_UMI_Simple --soloCBwhitelist Muraro.whitelist --soloCBstart 1 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen 4 --soloBarcodeReadLength 0 --soloCellFilter None --soloUMIfiltering MultiGeneUMI_CR --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --soloUMIdedup 1MM_CR --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCBmatchWLtype 1MM --readFilesCommand zcat

	# Cleanup
	mkdir /data/home/jgsm/VDS2/Mandrup/JM/Pancreas/Data/scRNAseq/Muraro/$z/
	mv Solo.out /data/home/jgsm/VDS2/Mandrup/JM/Pancreas/Data/scRNAseq/Muraro/$z/
	rm Aligned.out.sam
	mv Log* /data/home/jgsm/VDS2/Mandrup/JM/Pancreas/Data/scRNAseq/Muraro/$z/
	rm SJ.out.tab
	rm *.fq.gz*
	rm Donor
done
