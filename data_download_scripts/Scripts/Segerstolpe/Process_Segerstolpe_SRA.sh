# Load information about the number of samples in the study
Study="Segerstolpe.wget"
Donors=$(cut -f 1 $Study | sort | uniq)
Out="/work/57814/Data/scRNAseq/Segerstolpe"
Genome="/work/57814/Annotation/Processing/Package/hg38"

# Setup the files needed to be downloaded
cut -f 3 $Study | awk -F "/" '{ print $8 }' | paste $Study - | awk '{ print $1"\t"$2"\t"$5"\t"$4 }' > SRAdownload
nfiles=$(wc -l SRAdownload | tr ' ' '\t' | cut -f 1)
for (( c=1; c<=$nfiles; c++ ))
do
SRA=$( awk -v VAR=$c 'NR == VAR { print $3 }' SRAdownload)
FASTQ=$( awk -v VAR=$c 'NR == VAR { print $2 }' SRAdownload)
prefetch $SRA -C yes
fasterq-dump $SRA -o $FASTQ --include-technical -e 63
rm -rf $SRA
done

# Create a manifest file
VAR=$(ls *.gz)
for i in $VAR; do cell=$(echo $i | sed 's/.fastq//g'); echo -e $i"\t-\t"$cell >> manifest; done


# Run STAR
STAR --genomeDir $Genome --soloType SmartSeq --readFilesManifest ./manifest --soloUMIdedup Exact --soloStrand Unstranded --soloFeatures Gene GeneFull --outFilterScoreMin 30 --soloMultiMappers EM --soloCellFilter None --runThreadN 64 --outMultimapperOrder Random --outSAMmultNmax 1

# Cleanup
mv Solo.out $Out
rm Aligned.out.sam
mv Log* $Out
rm SJ.out.tab
rm *.fastq
rm Donor
rm SRAdownload
rm manifest
