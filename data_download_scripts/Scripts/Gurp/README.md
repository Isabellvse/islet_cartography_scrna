## Gurp

- `GEO`: GSE150724
- `BioProject`: PRJNA633387
- `PMID`: 35440614
- `DOI`: 10.1038/s41467-022-29588-8
- `Tissue`: Isolated islets
- `Year`: 2022

### Preprocessing scripts

- `process_Gurp.sh`: Download and align data with STARsolo
- `alignment_qc_Gurp.sh`: Combine all alignment qaulity control parameters for all samples in a table

### Dataset at a glance

- `n samples`: 3 Non diabetic
- `method`: Single Cell 3' v3
- `Library layout`: Paired

### STARsolo Alignment paramters

- `--soloType CB_UMI_Simple`: Sets the type of single-cell RNA-seq data, (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium.
- `--soloUMIlen 12`: UMI length is 12 bp
- `--soloBarcodeReadLength 0`: Length of barcode read. Not defined, do not check
- `--soloCellFilter None`: Disables cell filtering, meaning all barcodes will be considered.
- `--soloUMIfiltering MultiGeneUMI_CR`: Filters UMIs: basic + remove lower-count UMIs that map to more than one gene
- `--outMultimapperOrder Random`: n outputs multiple alignments for each read in random order, and also also randomizes the choice of the primary alignment from the highest scoring alignments.
- `--outSAMmultNmax 1`: Limits the number of alignments per read in the output SAM file to 1. Will output exactly one SAM line for each mapped read
- `--soloCBwhitelist Gurp.whitelist`: whitelist, the same as 3M-february-2018.txt
- `--soloUMIdedup 1MM_CR`: CellRanger2-4 algorithm for 1MM UMI collapsing.
- `--soloFeatures Gene GeneFull Velocyto`: Specifies the features to be counted, including genes (reads match the gene transcript), GeneFull (count all reads overlapping genes' exons and introns), and Velocyto (This option will calculate Spliced, Unspliced, and Ambiguous counts)
- `--outFilterScoreMin 30`: Sets the minimum score for filtering out low-quality reads.
- `--soloMultiMappers EM`: Uses the Maximum Likelihood Estimation for handling multi-mapping reads.
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`: Multiple matches in whitelist with 1 mismatched base or N mismatched bases allowed, posterior probability calculation is used choose one of the matches. Allowed cell barcodes have to have at least one read with exact match.
- `--clipAdapterType CellRanger4`: Uses the adapter clipping method similar to Cell Ranger version 4. 5' TSO (template switch oligo) adapter and 3' polyA-tail clipping of the reads to better match CellRanger 
