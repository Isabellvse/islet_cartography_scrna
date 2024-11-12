# Define variables
Genome="/data/home/jgsm/Pancreas/hg38/"
Out="/data/home/jgsm/VDS2/Mandrup/JM/Pancreas/HPAP/New/scRNAseq"
mkdir $Out
Study="PatchSeq.samples"
nfile=$(wc -l $Study | tr ' ' '\t' | cut -f 1)

# Loop
for ((i=1; i<=$nfile; i++ ))
do
	# Extract donor and chemistry information from the samples files (subset from official HPAP metadata)
	Donor=$(awk -v VAR=$i 'NR == VAR { print $1 }' $Study)

	# Download the data
	sftp hpapsftp@hpap-test.pmacs.upenn.edu:/hpapdata/"$Donor/Islet\ Studies/Islet\ physiology\ studies/Electrophysiology/Patch-Seq\ studies/Patch-Seq\ sequencing/"*fastq.gz .

	# Create a manifest file
	VAR=$(ls *.gz)
	for i in $VAR; do echo $i | sed "s/$Donor//g" - | sed 's/_patchseq_//g' - | sed 's/-R1_fastq-data.fastq.gz//g' - | sed 's/-R2_fastq-data.fastq.gz//g' - >> Cells; done
	cat Cells | sort | uniq -c | awk '$1 == "2" { print $2 }' - > ValidCells
	rm Cells
	VAR=$(cat ValidCells)
	for i in $VAR; do echo -e $Donor"_patchseq_"$i"-R1_fastq-data.fastq.gz\t"$Donor"_patchseq_"$i"-R2_fastq-data.fastq.gz\t"$i >> manifest; done
	rm ValidCells
	
	# Align the data using STARsolo
  	STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --soloCellFilter None --outFilterScoreMin 30 --soloMultiMappers EM --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

	# Move the results from STARsolo to VDS2
	mkdir $Out/$Donor
	mv Solo.out/ $Out/$Donor
	mv Log.* $Out/$Donor

	# Cleanup the processing folder
	rm *.fastq.gz
	rm manifest
	rm Aligned.out.sam
	rm SJ.out.tab
done
