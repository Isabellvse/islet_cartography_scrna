## Motakis

- `GEO`: GSE221156
- `BioProject`: PRJNA913127
- `PMID`: None
- `DOI`: None
- `Tissue`: Isolated islets
- `Year`: 2023

### Preprocessing scripts

- `process_Motakis.sh`: Download and align data with STARsolo
- `alignment_qc_Motakis.sh`: Combine all alignment qaulity control parameters for all samples in a table

### Dataset at a glance

- `n samples`: 17 Non diabetic, 14 Prediabetic and 17 type 2 diabetic
- `method`: Single-Cell RNA-seq, 10x Genomics Single Cell 3' v3
- `Library layout`: Paired

### STARsolo Alignment paramters

- `--soloType CB_UMI_Simple:` Sets the type of single-cell RNA-seq data, (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium.
- `--soloFeatures Gene GeneFull:` Specifies the features to be counted, including genes (reads match the gene transcript) and genefull (count all reads overlapping genes' exons and introns).
- `--soloCellFilter None:` Disables cell filtering, meaning all barcodes will be considered.
- `--soloUMIlen 12:` Sets the length of the UMI to 12 bases.
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts:` Multiple matches in whitelist with 1 mismatched base or N mismatched bases allowed, posterior probability calculation is used choose one of the matches. Allowed cell barcodes have to have at least one read with exact match. 
- `--clipAdapterType CellRanger4:` Uses the adapter clipping method similar to Cell Ranger version 4. 5' TSO (template switch oligo) adapter and 3' polyA-tail clipping of the reads to better match CellRanger 
- `--outFilterScoreMin 30:` Sets the minimum score for filtering out low-quality reads.
- `--soloMultiMappers EM:` Uses the Maximum Likelihood Estimation for handling multi-mapping reads.
- `--soloUMIfiltering MultiGeneUMI_CR:` Filters UMIs: basic + remove lower-count UMIs that map to more than one gene
- `--soloUMIdedup 1MM_CR:` CellRanger2-4 algorithm for 1MM UMI collapsing.
- `--outMultimapperOrder Random:` n outputs multiple alignments for each read in random order, and also also randomizes the choice of the primary alignment from the highest scoring alignments.
- `--outSAMmultNmax 1:` Limits the number of alignments per read in the output SAM file to 1. Will output exactly one SAM line for each mapped read
- `--soloCBwhitelist:` 3M-february-2018.txt
- `--soloBarcodeReadLength 0:` Length of barcode read. Not defined, do not check