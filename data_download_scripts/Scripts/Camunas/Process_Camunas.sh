# Conda environment
eval "$(conda shell.bash hook)"
conda activate /work/57814/conda/JM

# Define variables
Genome="/work/57814/Annotation/Processing/Package/hg38/"
Out="/work/57814/Data/Patchseq/Camunas/Preprocessed/"
Study="Camunas.wget"

# Download data
cut -f 4 $Study | xargs -n 1 -P 8 wget -q

# Create a manifest file
awk 'NR%2==1 { print $2"\t"substr($4, 58,12) }' Camunas.wget > SRRs
sed -i 's|_1||g' SRRs
sed -i 's|/||g' SRRs
awk '{ print $2"_1.fastq.gz\t"$2"_2.fastq.gz\t"$1 }' SRRs > manifest

# Run STAR
STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

# Cleanup
mv Solo.out $Out/
rm Aligned.out.sam
mv Log* $Out/
rm SJ.out.tab
rm *.fq.gz*
rm SRRs
rm manifest
