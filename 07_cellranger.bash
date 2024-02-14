#!/bin/bash
#
#SBATCH --job-name=cellranger
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --partition=basic
#SBATCH --output=log/cellr_count-%j.out
#SBATCH --error=log/cellr_count-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=june.ordonez@univie.ac.at
# #SBATCH --test-only

## SUMMARY: This script performs the initial part of analysing single-cell data: mapping the reads
## to the reference genome by cellranger

##load module
module load cellranger/7.1.0
module list

#go to the preferred working directory
WORK_DIR=/scratch/course/2023w300106/juny/day06
cd ${WORK_DIR}

#### Set the paths and variables (if any) ###
GENOME_DIR=/scratch/course/2023w300106/ #path to the directory containing the fastq files
GTF_DIR=/scratch/course/2023w300106/
OUT_DIR=/scratch/course/2023w300106/juny/day06
SCS_DIR=/scratch/course/2023w300106/12hr1

RUN_NAME=Nv2 ## Input run name
GENOME_FILE=Nv2_wnt4_pcna_fluo.fa
GTF_FILE=Nv2_wnt4_pcna_fluo.gtf
SCS_SAMPLE_NAME=89085

## ---------------------------------------------------------------------------------------------- ##
## Perform cellranger: the first part is the generation of indexed reference and the second part
## is the scRNAseq read mapping to the indexed reference
start_run=$(date +%M)

echo "Starting: cellranger.."

cellranger mkref --genome ${RUN_NAME} \
	  --fasta ${GENOME_FILE} \
	  --genes ${GTF+FILE} \
	  --nthreads 16

## ------------------------------------------------------------------------------------------------ ##
## IMPORTANT NOTES: If you have no idea what kind of chemistry was used to generate the scrnaseq data,
## leave the --chemistry flag because this automatically identifies it. But, for older-generated
## libraries, it is recommended to specify the chemistry (according to ACole)

cellranger count --id ${RUN_NAME}_count \
        --transcriptome ${OUT_DIR}/${RUN_NAME} \
        --fastqs ${SCS_DIR}/run1,${SCS_DIR}/run2,${SCS_DIR}/run3 \
        --sample ${SCS_SAMPLE_NAME} \
        --nosecondary \
        --localcores ${SLURM_CPUS_PER_TASK}

echo "Cellranger finished!"
end_run=$(date +%M)
echo "Elapsed time: $(($end_run-$start_run)) m"

echo " Check the ${RUN_NAME}_count/outs/web_summary.html for the summary of the output"