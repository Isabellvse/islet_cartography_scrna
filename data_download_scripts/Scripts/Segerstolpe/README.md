## Segerstolpe

- `Acession`: E-MTAB-5061
- `BioProject`: NA
- `PMID`: 27667667
- `DOI`: 10.1016/j.cmet.2016.08.020
- `Tissue`: FACS sorted islet cells
- `Year`: 2016

### Preprocessing scripts

- `process_Segerstolpe.sh`: Download and align data with STARsolo
- `alignment_qc_Segerstolpe.sh`: Combine all alignment qaulity control parameters for all samples in a table

### Dataset at a glance

- `n samples`:  6 Non diabetic,  4 type 2 diabetic
- `method`: Plate based, Smart-seq2
- `Library layout`: Paired

### STARsolo Alignment paramters

- `--soloType SmartSeq`: option produces cell/gene (and other features) count matrices, using rules similar to the droplet-based technologies. The differnces are (i) individual cells correspond to different FASTQ files,there are no Cell Barcode sequences, and "Cell IDs" have to be provided as input (ii) there are no UMI sequences, but reads can be deduplicated if they have identical start/end coordinates.
- `--readFilesManifest ./manifest`: A file manifest with a list of FASTQ fiels and cell ids
- `--soloUMIdedup Exact`: only exactly matching UMIs are collapsed
- `--soloStrand Unstranded`: only exactly matching UMIs are collapsed
- `--soloFeatures Gene GeneFull`: Specifies the features to be counted, including genes (reads match the gene transcript) and genefull (count all reads overlapping genes' exons and introns).
- `--outFilterScoreMin 30`: Sets the minimum score for filtering out low-quality reads.
- `--soloMultiMappers EM`: Uses the Maximum Likelihood Estimation for handling multi-mapping reads.
- `--soloCellFilter None`: do not output filtered cells
- `--outMultimapperOrder Random`: n outputs multiple alignments for each read in random order, and also also randomizes the choice of the primary alignment from the highest scoring alignments.
- `--outSAMmultNmax 1`: Limits the number of alignments per read in the output SAM file to 1. Will output exactly one SAM line for each mapped read
