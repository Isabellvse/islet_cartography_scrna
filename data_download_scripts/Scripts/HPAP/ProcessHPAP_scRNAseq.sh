# Conda environment
eval "$(conda shell.bash hook)"
conda activate /work/57814/conda/JM

# Define variables
Genome="/work/57814/Annotation/Processing/Package/hg38/"
Out="/work/57814/Data/scRNAseq/HPAP_10X/Preprocessed/"
Study="HPAP.samples_070223"
nfile=$(wc -l $Study | tr ' ' '\t' | cut -f 1)

# Loop
for ((i=1; i<=$nfile; i++ ))
do
	# Extract donor and chemistry information from the samples files (subset from official HPAP metadata)
	Donor=$(awk -v VAR=$i 'NR == VAR { print $1 }' $Study)
	Chemistry=$(awk -v VAR=$i 'NR == VAR { print $2 }' $Study)

	# Download the data
	sftp -i /work/57814/Annotation/Processing/Package/Scripts/HPAP_JM hpapsftp@hpap-test.pmacs.upenn.edu:/hpapdata/"$Donor/Islet\ Studies/Islet\ molecular\ phenotyping\ studies/Single-cell\ RNAseq/Upenn_scRNAseq/fastq/"*fastq.gz .

	# Cat the individual data files and remove the individual files
	cat *R1* > R1.fq.gz
	cat *R2* > R2.fq.gz
	rm *.fastq.gz
	
	# Align the data using STARsolo
	if [[ "$Chemistry" == "v2" ]];
	then
        	STAR --genomeDir $Genome --readFilesIn R2.fq.gz R1.fq.gz --soloType CB_UMI_Simple --soloFeatures Gene GeneFull Velocyto --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat --soloCBwhitelist HPAP.whitelist_v2 --soloBarcodeReadLength 0
	else
        	STAR --genomeDir $Genome --readFilesIn R2.fq.gz R1.fq.gz --soloType CB_UMI_Simple --soloFeatures Gene GeneFull Velocyto --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat --soloCBwhitelist HPAP.whitelist_v3 --soloBarcodeReadLength 0 --soloUMIlen 12
	fi

	# Move the results
	mkdir $Out/$Donor
	mv Solo.out/ $Out/$Donor
	mv Log.* $Out/$Donor

	# Cleanup the processing folder
	rm *.fq.gz
	rm Aligned.out.sam
	rm SJ.out.tab
done
