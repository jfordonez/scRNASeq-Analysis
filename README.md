# single-cell RNAseq analysis: a pre-analysis pipeline
The repository contains scripts used in single-cell analysis course (Jan 8 - Feb 5, 2023)

To summarize, the scripts provide a step-by-step commands on how to build a reference gene-model dataset that will help guide the mapping of scRNA-seq reads. Remember, we are more interested in reads that mapped to genes because we are concerned with the differential gene expression profiles of cells.

## What these scripts do:
1. Apply quality control (trim and filter) to bulk RNAseq reads and then map them to a reference genome (scripts 01 and 02)
2. Predict genes and create an annotation file (gff/gtf) that contains the position and configuration of each gene along the genome based on the RNAseq mapping. (scripts 03 and 04)
3. Annotate the predicted genes: what are they and what their putative function based on conserved orthologous protein domains. (script 06)
4. Preliminary analysis and quality check of the scRNAseq data by mapping them to the genome that features the guide gtf generated from scripts 01 to 06. (script 07)

## What these scripts don't do:
- post-analyze scRNAseq datasets (no UMAP o tSNE)


## Data used in the example dataset (accession nos.:
1. Nematostella vectensis RNAseq reads (blastula stage): SRR24348406
2. Nematostella vectensis Genome: GCF_932526225
3. For cellranger = gtf and genome files were provided by the instructors


