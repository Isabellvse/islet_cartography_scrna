# Load information about the number of samples in the study
Study="Segerstolpe.wget"
Donors=$(cut -f 1 $Study | sort | uniq)
Genome="/work/57814/Annotation/Processing/Package/hg38"
Out="/work/57814/Data/scRNAseq/Segerstolpe"
mkdir -p $Out

# Download all files
cut -f 4 $Study | xargs -n 1 -P 8 wget -q

# Create a manifest file
awk '{ print substr($4, 58, 60)"\t-\t"$2 }' Segerstolpe.wget > manifest

# Run STAR
STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 64 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat

# Cleanup
mv Solo.out $Out/
rm Aligned.out.sam
mv Log* $Out/
rm SJ.out.tab
rm *.fq.gz*
rm manifest
