# Load information about the number of samples in the study
Study="Muraro.wget"
Donors=$(cut -f 1 $Study | sort | uniq)
Out="/work/57814/Data/scRNAseq/Muraro"
Genome="/work/57814/Annotation/Processing/Package/hg38"

# Setup the files needed to be downloaded
cut -f 3 $Study | awk -F "/" '{ print $8 }' | paste $Study - | awk '{ print $1"\t"$2"\t"$5"\t"$4 }' > SRAdownload
nfiles=$(wc -l SRAdownload | tr ' ' '\t' | cut -f 1)
for (( c=1; c<=$nfiles; c++ ))
do
SRA=$( awk -v VAR=$c 'NR == VAR { print $3 }' SRAdownload)
DONOR=$( awk -v VAR=$c 'NR == VAR { print $1 }' SRAdownload)
sratoolkit.3.0.2-ubuntu64/bin/prefetch $SRA -C yes --output-file $SRA.srr
sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump $SRA.srr --include-technical -e 63
rm -rf $SRA.srr

# Run STAR
STAR --genomeDir $Genome --readFilesIn *_2.fastq *_1.fastq --soloType CB_UMI_Simple --soloCBwhitelist Muraro.whitelist --soloCBstart 1 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen 4 --soloBarcodeReadLength 0 --soloCellFilter None --soloUMIfiltering - --runThreadN 20 --outMultimapperOrder Random --outSAMmultNmax 1 --soloUMIdedup Exact --soloFeatures Gene GeneFull Velocyto --outFilterScoreMin 30 --soloMultiMappers EM --soloCBmatchWLtype 1MM

# Cleanup
mkdir $Out/Preprocessed/$DONOR
mv Solo.out $Out/Preprocessed/$DONOR/
rm Aligned.out.sam
mv Log* $Out/Preprocessed/$DONOR
rm SJ.out.tab
rm *.fastq
done
rm SRAdownload
