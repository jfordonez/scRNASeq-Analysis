#!/bin/bash
#
#SBATCH --job-name=Trim&Map
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --partition=basic
#SBATCH --output=log/trimmo-%j.out
#SBATCH --error=log/trimmo-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=mail.ac.at
# #SBATCH --test-only

#Load modules needed
module load conda/miniconda3
module load trimmomatic/0.39
module load star/2.7.11a
module list

#go to the preferred working directory
WORK_DIR=/scratch/course/2023w300106/juny/day02/01_Trimmomatic
cd ${WORK_DIR}

### Set the variables (if any) ###
RUN_NAME=Evec_SRR24348406
READ1=SRR24348406_1
READ2=SRR24348406_2

## Set the paths of files to be used
GENOME_DIR=/scratch/course/2023w300106/juny/day01/01_dataset/ncbi_dataset/data/GCF_932526225.1/ #use full path
READ_DIR=/scratch/course/2023w300106/juny/day01/01_dataset/SRR24348406  #path to the directory containing the fastq files
OUT_DIR=/scratch/course/2023w300106/juny/day02/01_Trimmomatic
ILL_CLIP=/scratch/course/2023w300106/juny/day02/01_Trimmomatic/adapters.fa ## guide ref seq for illumina adapters

## --------------------------------------------------------------------------------------- ##
### 1. Run quality control using trimmomatic ###

echo "Starting QC ...."
echo "Starting: Trimmomatic"
start_run=$(date +%M) # start run timer

trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} -trimlog trimmomatic_${RUN_NAME}.log \
               -summary trim.summary.${RUN_NAME}.txt \
               ${READ_DIR}/${READ1}.fastq ${READ_DIR}/${READ2}.fastq \
               ${READ1}.paired.fq ${READ1}.unpaired.fq \
               ${READ2}.paired.fq ${READ2}.unpaired.fq \
               ILLUMINACLIP:${ILL_CLIP}:2:30:10 \
               SLIDINGWINDOW:6:15 MINLEN:75

end_run=$(date +%M)
echo "Trimmomatic finished"
echo "Elapsed time: $(($end_run-$start_run)) m"

## ----------------------------------------------------------------------------------------- ##
### 2. Run alignment of trimmed files to the reference genome ###
start_run=$(date +%M)

echo "Starting: alignment of trimmed reads with STAR..."

STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
        --genomeDir ${GENOME_DIR} \
        --readFilesIn ${OUT_DIR}/${READ1}.paired.fq ${OUT_DIR}/${READ2}.paired.fq \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outFileNamePrefix ${RUN_NAME}_trim_


end_run=$(date +%M)
echo "STAR finished"
echo "Elapsed time: $(($end_run-$start_run)) m"

conda deacivate
ml purge