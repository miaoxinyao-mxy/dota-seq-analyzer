# ARGMapper

ARGMapper identifies the antibiotic resistance genes \(ARGs\) present  
in each species of bacteria from a sample.

This software is intended to be used for bacteria samples of unknown  
composition. Single-cell sequencing with short, paired reads should  
be used. Further, this software was designed for use with the rest  
of the DoTA-Seq workflow.

Primary use case: by running this software on the same bacterial sample  
at two different timepoints, one can observe the horizontal gene  
transfer \(HGT\).

## Quick Start

### System Requirements

### Installation & Setup

### Command

Required Inputs:
\(provided you're using the default of "yes"  
to both ASVs and sub-ARGs\)
- R1 FASTQ
- R2 FASTQ
- Primers CSV
- Sub-ARG database FASTA  
\(to be used for BLAST; omit if not using sub-ARGs\)
- Output folder

Run the following command:
```
argmapper --r1_fastq example_r1.fastq --r2_fastq example_r2.fastq  
--primers example_primers.csv --db_sub_arg_fasta example_db.fasta  
-o output_folder
```

The software will print out the paths for the final output files.

### Interpreting Output

#### ASV-ARG Heat Map

Each square shows the % of cells identified as belonging to the  
given ASV, that have the given sub-ARG.

#### Barcode Summary TSV

This shows all MLE taxonomy, ASV, & sub-ARG information on a  
per-cell basis \(note that each barcode is assumed to be associated  
with only one cell, and vice versa\)

## Software Description

Here we describe the general stages of ARGMapper's pipeline.

First, Kraken2 is used for initial 16s taxonomic classifications.  
The software then organizes and processes the data into a  
barcode-based summary. Next, ASV \(Amplicon Sequence Variant\)  
information, used for taxonomic classification, is appended to this  
summary. The software also splits the ARG columns from the barcode summary  
into sub-ARGs, and then uses BLAST to identify names for each sub-ARG sequence.

Statistical algorithms are used throughout the pipeline, for  
filtering and classification purposes. The final ASV-ARG heat map  
output shows the magnitude of each sub-ARG's presence in each ASV  
taxonomic classification of bacteria.

For more information on how ARGMapper works, please reference our paper.

## Optional Parameters

Optional parameters for ARGMapper include:
1) `--threads` 
2) `--use_asvs`
3) `--use_sub_args`
4) `--filter_taxonomic_classifications`
5) `--filter_low_confidence_single_asvs`

There are also other parameters that don't require changes as often.  
These can be modified directly in the "driver.sh" file.
These parameters include maximum shifts and/or mismatch values  
\(for comparing primers, barcodes, and sub-ARG sequences\),  
filtering limits, and statistical parameters.

## Citation

## Contact Information