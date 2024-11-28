# Conda environment
eval "$(conda shell.bash hook)"
conda activate /work/57814/conda/JM

# Define variables
folders=$(ls -d SRR18593*)

study_name="Wang_Sander"
Study="${study_name}.wget"
Out="/work/scRNAseq/${study_name}/Preprocessed"
mkdir -p "$Out"
Donors=$(cut -f 1 "$Study" | sort | uniq)
Genome="/work/islet_cartography_scrna/data_download_scripts/hg38/"
whitelist="/work/islet_cartography_scrna/whitelist/737K-arc-v1.txt.txt"

for z in $Donors; do
    # Setup the files needed to be downloaded for each donor
    awk -v VAR=$z '$1 == VAR { print $0 }' "$Study" > Donor
    awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' Donor > Download

    # Download using aria2
    aria2c -i Download -j10 --check-integrity=true --save-session failed.downloads
    has_error=$(wc -l < failed.downloads)
    
    while [ $has_error -gt 0 ]; do
        echo "Still has $has_error errors, rerun aria2 to download ..."
        mv failed.downloads Download
        aria2c -i Download -j10 --check-integrity=true -c --save-session failed.downloads
        has_error=$(wc -l < failed.downloads)
        sleep 10
    done

    # Run STAR 
            STAR --genomeDir $Genome --readFilesIn $i/*R2*.fastq.gz $i/*R1*.fastq.gz --soloType  \
            CB_UMI_Simple --soloFeatures Gene GeneFull Velocyto --soloCellFilter None --soloCBmatchWLtype \
            1MM_multi_Nbase_pseudocounts --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloMultiMappers EM \
            --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --runThreadN 60 --outMultimapperOrder Random \
            --outSAMmultNmax 1 --readFilesCommand zcat --soloCBwhitelist $whitelist --soloBarcodeReadLength 0 --soloUMIlen 12

    # Cleanup: Move results to donor-specific folder
    rm failed.downloads
	mkdir $Out/$z/
	mv Solo.out $Out/$z/
	rm Aligned.out.sam
	mv Log* $Out/$z/
	rm SJ.out.tab
	rm *.fastq
	rm Donor
	rm Download
	rm *.sra
done
