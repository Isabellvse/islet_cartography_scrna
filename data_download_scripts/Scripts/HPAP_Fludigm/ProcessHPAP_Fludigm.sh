# Define variables
Genome="/work/57814/Annotation/Processing/Package/hg38"
Out="/work/57814/Data/scRNAseq/HPAP_Fluidigm"
mkdir $Out
Study="Fludigm.samples"
nfile=$(wc -l $Study | tr ' ' '\t' | cut -f 1)

# Loop
for ((ids=1; ids<=$nfile; ids++ ))
do
	# Extract donor and chemistry information from the samples files (subset from official HPAP metadata)
	Donor=$(awk -v VAR=$ids 'NR == VAR { print $1 }' $Study)

	# Download the data
	sftp -i /home/ucloud/.ssh/id_rsa -r hpapsftp@hpap-test.pmacs.upenn.edu:/hpapdata/"$Donor/Islet\ Studies/Islet\ molecular\ phenotyping\ studies/Single-cell\ RNAseq/Upenn_scRNAseq/fastq/"*_scRNA*gz .

	# Check for errors
	Files=$(ls *.gz)
	touch Err.file
	for i in $Files; do gzip -t $i 2>&1 | awk 'NR > 1 { print $0 }' - | awk -F ":" '{ print $2 }' - | sed 's/ //g' >> Err.file; done
	Del=$(cat Err.file)
	for i in $Del; do rm $i; done

	# Create a manifest file
	VAR=$(ls *.gz)
	for i in $VAR; do echo $i | sed "s/$Donor//g" | sed 's/_scRNA_//g' | sed 's/_fastq-data.fastq.gz//g' >> Cells; done
	VAR=$(cat Cells)
	for i in $VAR; do echo -e $Donor"_scRNA_"$i"_fastq-data.fastq.gz\t-\t"$i >> manifest; done
	rm Cells
	
	# Align the data using STARsolo
	STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --soloCellFilter None --outFilterScoreMin 30 --soloMultiMappers EM --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

	# Move the results from STARsolo to VDS2
	mkdir $Out/$Donor
	mv Solo.out/ $Out/$Donor
	mv Log.* $Out/$Donor
	mv Err.file $Out/$Donor

	# Cleanup the processing folder
	rm *.gz
	rm Aligned.out.sam
	rm SJ.out.tab
	rm manifest
done
