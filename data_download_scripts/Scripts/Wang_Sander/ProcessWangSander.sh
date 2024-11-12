# Conda environment
eval "$(conda shell.bash hook)"
conda activate /work/57814/conda/JM

# Define variables
Genome="/work/57814/Annotation/Processing/Package/hg38/"
Out="/work/57814/Data/scRNAseq/Wang_Sander/Preprocessed/"
folders=$(ls -d SRR18593*)

# Loop
for i in $folders
do
	# Define the donor ID
	Donor=$(ls $i | tr '_' '\t' | awk 'NR == 1 { print $1"_"$2 }')

	# Align the data using STARsolo
        STAR --genomeDir $Genome --readFilesIn $i/*R2*.fastq.gz $i/*R1*.fastq.gz --soloType CB_UMI_Simple --soloFeatures Gene GeneFull Velocyto --soloCellFilter None --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 60 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat --soloCBwhitelist 737K-arc-v1.txt --soloBarcodeReadLength 0 --soloUMIlen 12

	# Move the results
	mkdir $Out/$Donor
	mv Solo.out/ $Out/$Donor
	mv Log.* $Out/$Donor

	# Cleanup the processing folder
	rm Aligned.out.sam
	rm SJ.out.tab
done
