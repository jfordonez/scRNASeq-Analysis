#!/bin/bash
#
#SBATCH --job-name=reads2transciptome
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --partition=basic
#SBATCH --output=log/STAR_trans-%j.out
#SBATCH --error=log/STAR_trans-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=june.ordonez@univie.ac.at
# #SBATCH --test-only

## SUMMARY: Align the RNAseq reads to a draft transcriptome available in ncbi and assess
## the differences between mapping to a genome.

##load modules
module load star/2.7.11a
module list

#go to the preferred working directory
WORK_DIR=/scratch/course/2023w300106/juny/day04/02_Nvec_Transcriptome
cd ${WORK_DIR}

## Set the paths and variables (if any) ##
READ_DIR=/scratch/course/2023w300106/juny/day02/01_Trimmomatic/  #path to the directory containing the fastq files
TRANSC_DIR=/scratch/course/2023w300106/juny/day04/02_Nvec_Transcriptome/HADO01
OUT_DIR=/scratch/course/2023w300106/juny/day04/02_Nvec_Transcriptome

RUN_NAME=Evec_SRR24348406 ##Input run name
READ1=SRR24348406_1
READ2=SRR24348406_2
TRANSC_FILE=HADO01.fasta #This was downloaded from NCBI database

### Run STAR . . .#
start_run=$(date +%M)

echo "Starting: STAR mapping. . ."

STAR --runMode=genomeGenerate \
      --runThreadN ${SLURM_CPUS_PER_TASK} \
      --genomeDir ${TRANSC_DIR} \
      --genomeFastaFiles ${TRANSC_DIR}/${TRANSC_FILE} \
      --genomeSAindexNbases 11

STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
        --genomeDir ${TRANSC_DIR} \
        --readFilesIn ${READ_DIR}/${READ1}.paired.fq ${READ_DIR}/${READ2}.paired.fq \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outFileNamePrefix ${RUN_NAME}_trans_


echo "STAR mapping ended.."


end_run=$(date +%M)
echo "Elapsed time: $(($end_run-$start_run)) m"
