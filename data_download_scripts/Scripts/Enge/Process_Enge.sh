# Load information about the number of samples in the study
Study="Enge.wget"
Genome="/work/57814/Annotation/Processing/Package/hg38"
Donors=$(cut -f 1 $Study | sort | uniq)
Out="/work/57814/Data/scRNAseq/Enge/Preprocessed"
mkdir $Out

# Loop over donors
for z in $Donors; do
	# Setup the files needed to be downloaded
	awk -v VAR=$z '$1 == VAR { print $0 }' $Study > Donor		
	awk '{ print $4"\n out="$3"\n checksum=md5="$5 }' Donor > Download

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

	# Create the manifest file
	VAR=$(ls *.gz)
	for i in $VAR; do echo $i | sed "s/_R1.fq.gz//g" - | sed "s/_R2.fq.gz//g" >> Cells; done
	cat Cells | sort | uniq -c | awk '$1 == "2" { print $2 }' - > ValidCells
	rm Cells
	VAR=$(cat ValidCells)
	for i in $VAR; do echo -e $i"_R1.fq.gz\t"$i"_R2.fq.gz\t"$i >> manifest; done
	rm ValidCells

	# Run STAR
	STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

	# Cleanup
	mkdir $Out/$z/
	mv Solo.out $Out/$z/
	rm Aligned.out.sam
	mv Log* $Out/$z/
	rm SJ.out.tab
	rm *.fq.gz*
	rm Donor
	rm Download
	rm manifest
	mv failed.downloads $Out/$z/
done
