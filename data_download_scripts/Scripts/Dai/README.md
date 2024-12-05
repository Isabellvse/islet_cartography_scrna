## Dai

- `GEO`: GSE164875
- `BioProject`: PRJNA692244
- `PMID`: 35108513
- `DOI`: 10.1016/j.cmet.2021.12.021
- `Tissue`: Pancreas cells
- `Year`: 2021

### Preprocessing scripts

- `process_Dai.sh`: Download and align data with STARsolo
- `alignment_qc_Dai.sh`: Combine all alignment qaulity control parameters for all samples in a table

### Dataset at a glance

- `n samples`:  640 pancreatic cells from Non diabetic and type 2 diabetic individuals
- `method`: Patch-seq data, aligned with Smart-seq2
- `Library layout`: paired-end reads (100 bp)

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

#### Example of fastq file:

Read 1:

@SRR13439816.1 1/1
GTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF::FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:

Read 2: 

@SRR13439816.1 1/2
GCTTCCCATGTACTCTGCGTTGATACCACTGCTTCCCATGTACTCTGCGTTGATACCACTGCTTCCCATGTACTCTGCGTTGATACCACTGCTTCCCATG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FF:FFFFFFFFFFF:FFFFFFFFF:F,,,:FF:FFFFFFFFFFF,FFFFF:FF