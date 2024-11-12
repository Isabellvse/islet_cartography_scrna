# Load information about the number of samples in the study
Study="Enge.wget"
Genome="/work/57814/Annotation/Processing/Package/hg38"
Donors=$(cut -f 1 $Study | sort | uniq)
Out="/work/57814/Data/scRNAseq/Enge/Preprocessed"
mkdir $Out

# Create SRAfiles
awk 'NR % 2 == 1 { print $0 }' $Study | cut -f 4 | awk -F "/" '{ print $8 }' > SRAids
awk 'NR % 2 == 1 { print $0 }' $Study | paste - SRAids | awk '{ print $1"\t"$2"\t"$6 }' > SRAfiles

# Loop over donors
for z in $Donors; do
	# Setup the files needed to be downloaded
	awk -v VAR=$z '$1 == VAR { print $0 }' SRAfiles > SRAdownload	

	# Download using prefect and dump using fasterq-dump
	nfiles=$(wc -l SRAdownload | tr ' ' '\t' | cut -f 1)
	for (( c=1; c<=$nfiles; c++ ))
	do
		SRA=$( awk -v VAR=$c 'NR == VAR { print $3 }' SRAdownload)
		FASTQ=$( awk -v VAR=$c 'NR == VAR { print $2 }' SRAdownload)
		prefetch $SRA -C yes
		fasterq-dump $SRA -o $FASTQ --include-technical -e 63
		rm -rf $SRA
	done

	# Create the manifest file
	VAR=$(ls *.fastq)
	for i in $VAR; do echo $i | sed "s/_1.fastq//g" - | sed "s/_2.fastq//g" >> Cells; done
	cat Cells | sort | uniq -c | awk '$1 == "2" { print $2 }' - > ValidCells
	rm Cells
	VAR=$(cat ValidCells)
	for i in $VAR; do echo -e $i"_1.fastq\t"$i"_2.fastq\t"$i >> manifest; done
	rm ValidCells

	# Run STAR
	STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 64 --outMultimapperOrder Random --outSAMmultNmax 1

	# Cleanup
	mkdir $Out/$z/
	mv Solo.out $Out/$z/
	rm Aligned.out.sam
	mv Log* $Out/$z/
	rm SJ.out.tab
	rm *.fastq
	rm Donor
	rm SRAdownload
	rm manifest
done

# Final cleanup
rm SRAfiles
rm SRAids
