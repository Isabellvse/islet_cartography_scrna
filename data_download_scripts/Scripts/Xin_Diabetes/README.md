## Xin_Diabetes

- `GEO`: GSE114297
- `BioProject`:	PRJNA470834
- `PMID`: 29950394, 30380031, 33529174
- `DOI`: 10.2337/db18-0365, 10.1210/en.2018-00833, 10.1172/jci.insight.141553
- `Tissue`: Isolated islets
- `Year`: 2018

### Preprocessing scripts

- `process_Xin_Diabetes.sh`: Download and align data with STARsolo
- `alignment_qc_Xin_Diabetes.sh`: Combine all alignment qaulity control parameters for all samples in a table

### Dataset at a glance

- `n samples`: 12 Non diabetic
- `method`: Single-Cell RNA-seq, 10x Genomics Single Cell 3' v1
- `Library layout`: Paired

### STARsolo Alignment paramters
  
- `--soloType CB_UMI_Simple`: Sets the type of single-cell RNA-seq data, (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium.
- `--soloCBstart 1`: Cell barcode start base at bp 1
- `--soloCBlen 14`: Cell barcode length is 14 bp
- `--soloUMIstart 15`: UMI start is bp 13
- `--soloUMIlen 10`: UMI length is 10 bp
- `--soloBarcodeReadLength 1`: barcode read length should be equal to sum of soloCBlen+soloUMIlen
- `--soloCellFilter None`: Disables cell filtering, meaning all barcodes will be considered.
- `--soloUMIfiltering MultiGeneUMI_CR`: Filters UMIs: basic + remove lower-count UMIs that map to more than one gene
- `--outMultimapperOrder Random`: n outputs multiple alignments for each read in random order, and also also randomizes the choice of the primary alignment from the highest scoring alignments.
- `--outSAMmultNmax 1`: Limits the number of alignments per read in the output SAM file to 1. Will output exactly one SAM line for each mapped read
- `--soloCBwhitelist Xin.whitelist`: The whitelist, the same as 737K-april-2014_rc.txt
- `--soloUMIdedup 1MM_CR`: CellRanger2-4 algorithm for 1MM UMI collapsing.
- `--soloFeatures Gene GeneFull Velocyto`: Specifies the features to be counted, including genes (reads match the gene transcript), GeneFull (count all reads overlapping genes' exons and introns), and Velocyto (This option will calculate Spliced, Unspliced, and Ambiguous counts)
- `--outFilterScoreMin 30`: Sets the minimum score for filtering out low-quality reads.
- `--soloMultiMappers EM`: Uses the Maximum Likelihood Estimation for handling multi-mapping reads.
- `--soloCBmatchWLtype 1MM_multi_Nbase`: Filters UMIs: basic + remove lower-count UMIs that map to more than one gene