#!/bin/bash
#
#SBATCH --job-name=STAR
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --partition=basic
#SBATCH --output=log/mapping-%j.out
#SBATCH --error=log/mapping-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=mail@univie.ac.at
# #SBATCH --test-only

# Load modules needed
module load conda/miniconda3
module load star/2.7.11a
module list

# go to the preferred working directory
cd /scratch/course/2023w300106/juny/day01/02_Mapping

### Set the variables (if any) ###
RUN_NAME=Evec_blastula_SRR24348406
READ1=SRR24348406_1 #fastq read1
READ2=SRR24348406_2 #fastq read2

## Set the paths of files to be used. Use full path. ##
GENOME_DIR=GCF_932526225.1     
READ_DIR=/scratch/course/2023w300106/juny/day01/01_dataset/SRR24348406  #path to the directory containing the fastq files
OUR_DIR=/scratch/course/2023w300106/juny/day01/02_Mapping

## -------------------------------------------------------------------------------------------- ##
### 1. Generate index for the genome of E. vectensis with STAR ###
start_run=$(date +%M) #start run timer

## IMPORTANT NOTE: --genomeSAindexNbases set to 13 based on the genome size. Check manual for
## size guide.
echo "Generating reference index..."

STAR --runMode genomeGenerate \
       --runThreadN 8 \
       --genomeDir Evect_Star_indx \
       --genomeFastaFiles GCF_932526225.1/GCF_932526225.1_jaNemVect1.1_genomic.fna \
       --genomeSAindexNbases 13

## --------------------------------------------------------------------------------------------- ##
### 2. Run the STAR alignment ###

## IMPORTANT NOTE: --twopassMode Basic allows for higher read mapping to the novel junctions detected
## see manual for more information.

echo "Starting STAR alignment ...."

STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
        --genomeDir ${GENOME_DIR} \
        --readFilesIn ${READ_DIR}/${READ1}.fastq ${READ_DIR}/${READ2}.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outFileNamePrefix ${RUN_NAME}

end_run=$(date +%M)
echo "STAR mapping finished!!"
echo "Elapsed time: $(($end_run-$start_run)) m"

conda deactivate
ml purge