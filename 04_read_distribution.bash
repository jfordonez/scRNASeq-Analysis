#!/bin/bash
#
#SBATCH --job-name=featCnts
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=4G
#SBATCH --partition=basic
#SBATCH --output=log/featurecounts-%j.out
#SBATCH --error=log/feautrecounts-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=june.ordonez@univie.ac.at
# #SBATCH --test-only

## SUMMARY: This script will use featurecounts to get a handle of the read mapping of the RNAseq reads
## distribution on the genome

## load modules
module load subread/2.0.6
module list

## go to the preferred working directory
WORK_DIR=/scratch/course/2023w300106/juny/day04/01_TPM
cd ${WORK_DIR}

##  Set the paths for the files to be used (if any) ###
RUN_NAME=Evec_SRR24348406 ## Input run name
BAM_FILE=/scratch/course/2023w300106/juny/day03_GenomeAnnotation/01_StringTie/SRR24348406.f.bam ## the BAM file from star
GTF_FILE=/scratch/course/2023w300106/jmontenegro/ex2/annotation/tmp.gtf

## 1. Perform featurecounts 
start_run=$(date +%M)

echo "Starting: featureCounts . . ."

## IMPORTANT NOTES: 
## -C specifies that multimappers are not included in the analysis
## -p specifies that the libaries contain paired-end reads
## -s specifies that the libraries are strand specific

featureCounts -C -p -s 2 -Q 20 -T ${SLURM_CPUS_PER_TASK} \
        --countReadPairs \
        -a ${GTF_FILE} \
        -o ${RUN_NAME}.counts.tsv \
        ${BAM_FILE}

echo "featurecounts ended.."


end_run=$(date +%M)
echo "Elapsed time: $(($end_run-$start_run)) m"