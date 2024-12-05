## Camunas

- `GEO`: GSE124742 
- `BioProject`: PRJNA513339
- `PMID`: 32302527
- `DOI`: 10.1016/j.cmet.2020.04.005
- `Tissue`: Islet cells
- `Year`: 2020

### Preprocessing scripts

- `process_Camunas.sh`: Download and align data with STARsolo
- `alignment_qc_Camunas.sh`: Combine all alignment qaulity control parameters for all samples in a table

### Dataset at a glance

- `n samples`:  Non diabetic,  type 2 diabetic
- `method`: Patch-seq data, aligned with Smart-seq2
- `Library layout`: Paired

### STARsolo Alignment paramters

- `--soloType CB_UMI_Simple:` Sets the type of single-cell RNA-seq data, (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium.
- `--soloCBstart 1:` Cell barcode start base at position 1
- `--soloCBlen 12:` Cell barcode length is 12 bases
- `--soloUMIstart 13:` UMI start base is base 13
- `--soloUMIlen 8:` UMI length is 8 bases
- `--soloBarcodeReadLength 0:` Do not check length of the barcode read
- `--soloCellFilter None:` Disables cell filtering, meaning all barcodes will be considered.
- `--soloUMIfiltering MultiGeneUMI_CR:` Filters UMIs: basic + remove lower-count UMIs that map to more than one gene
- `--soloFeatures Gene GeneFull Velocyto:` Specifies the features to be counted, including genes (reads match the gene transcript) and genefull (count all reads overlapping genes' exons and introns).
- `--outMultimapperOrder Random:` n outputs multiple alignments for each read in random order, and also also randomizes the choice of the primary alignment from the highest scoring alignments.
- `--outSAMmultNmax 1:` Limits the number of alignments per read in the output SAM file to 1. Will output exactly one SAM line for each mapped read
- `--soloCBwhitelist none:` no whitelist
- `--outFilterScoreMin 30:` Sets the minimum score for filtering out low-quality reads.
- `--soloMultiMappers EM:` Uses the Maximum Likelihood Estimation for handling multi-mapping reads.
- `--soloUMIdedup 1MM_CR:` CellRanger2-4 algorithm for 1MM UMI collapsing.