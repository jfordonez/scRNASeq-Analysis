# single-cell RNAseq analysis: a pre-analysis pipeline

The repository contains scripts used in a single-cell analysis course taight by J. Montenegro and A. Cole in Univ. of Vienna (January 8 to February 5, 2023).

In summary, these batch scripts provide a step-by-step process for building a reference gene-model dataset and then mapping the scRNA-seq reads to a reference genome featuring the gene-model information initially generated. It is important to note that we are primarily interested in reads that map to the protein-coding regions of the genome, as we are primarily concerned with the differential gene expression profiles of cells.

## What these scripts do:
1. Apply quality control (trimming and filtering) to bulk RNAseq reads and then map them to a reference genome using scripts 01 and 02.
2. Predict genes and create an annotation file (gff/gtf) that contains the position and configuration of each gene along the genome based on the RNAseq mapping using scripts 03 and 04.
3. Annotate the predicted genes, including their putative function based on conserved orthologous protein domains, using script 06.
4. Perform a preliminary analysis and quality check of the scRNAseq data by mapping them to the genome that features the guide gtf generated from scripts 01 to 06 using script 07.

## What these scripts don't do:
- Perform post-analysis of scRNAseq datasets, such as UMAP or tSNE.

## Data used in the example dataset (ncbi accession numbers):
1. Nematostella vectensis RNAseq reads (blastula stage): SRR24348406
2. Nematostella vectensis Genome: GCF_932526225
3. For cellranger: the gtf and genome files were provided by the instructors.

