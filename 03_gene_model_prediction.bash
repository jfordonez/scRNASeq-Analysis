#!/bin/bash
#
#SBATCH --job-name=Annot
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=basic
#SBATCH --output=log/annotation-%j.out
#SBATCH --error=log/annotation-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=mail@univie.ac.at
# #SBATCH --test-only


## SUMMARY: This script works on the alignment files (.bam) and will be used to generate
## gene models, and ultimately a gtf file. Here, two approaches can be used.

#Load modules needed
module load conda/miniconda3
module load stringtie/2.2.1
module load samtools
module list

#go to the preferred working directory
cd /scratch/course/2023w300106/juny/day03_GenomeAnnotation/01_StringTie

#### Set the variables (if any) ###
RUN_NAME=Evec_SRR24348406 ## Input run name
BAM_FILE=Evec_SRR24348406_trim_Aligned.sortedByCoord.out.bam ## the BAM file out from star
FILTERED_BAM_FILE=SRR24348406_sort_filtr.bam ## sorted and filtered BAM file
GENOME_FILE=GCF_932526225.1_jaNemVect1.1_genomic.fna ##the genome file (in .fa/.fna)

### Set the paths of files to be used. Recommended: use full path
GENOME_DIR=/scratch/course/2023w300106/juny/day01/01_dataset/ncbi_dataset/data/GCF_932526225.1        
IN_DIR=/scratch/course/2023w300106/juny/day02/02_STAR_align #dir of your bam files
OUT_DIR_STRINGTIE=/scratch/course/2023w300106/juny/day03_GenomeAnnotation/01_StringTie #dir for the output

## ----------------------------------------------------------------------------------------------------- ##
echo "Starting gene model prediction ...."

## 1. Preprocess bam file to only contain high-quality mapping reads.

### IMPORTANT NOTES: filter out multimapping reads (0x4) and retain reads with mapQ = 20
### and proper-pair (0x2)

echo "Starting: BAM Filtering.........(^_______^ y)"
start_run=$(date +%M)

samtools view --bam \
      --min-MQ 20 \
      --require-flags 0x2 \
      --exclude-flags 0x4 \
      --threads ${SLURM_CPUS_PER_TASK} \
      --output ${OUT_DIR}/${FILTERED_BAM_FILE} \
      ${IN_DIR}/${BAM_FILE}

end_run=$(date +%M)
echo "BAM filtering finished"
echo "Elapsed time: $(($end_run-$start_run)) m"

## ----------------------------------------------------------------------------------------------------- ##
## The following sections are the tools to generate the general transfer format (.gtf), which provides
## some information on the predicted spatial configuration of genes in the genome. This also includes 
## position along the chromosome, exon-intron junctions, UTRs, intergentic regions
## etc. This is one of the core datasets that is essential for most downstream analyses

### 2. Run StingTie to predict gene models ###
start_run=$(date +%M)

echo "Starting: StringTie..."
stringtie ${OUT_DIR}/${FILTERED_BAM_FILE} \
      -l ${RUN_NAME}_StringTie \
      -o ${OUT_DIR_STRINGTIE}/${RUN_NAME}_StringTie.gtf \
      -p ${SLURM_CPUS_PER_TASK}
	  
end_run=$(date +%M)
echo "StringTie finished"
echo "Elapsed time: $(($end_run-$start_run)) m"

## ------------------------------------------------------------------------------------------------------ ##
### 3. Run Augustus to predict gene models ###

module load conda/miniconda3
conda activate augustus-3.5.0
echo "Conda env activated: $CONDA_DEFAULT_ENV"

## SET NAMES ##
CHROM_NAME=NC_06044.1 ## choose any chromosome of interest. This should be in ".fasta" name ending

## update working dir ##
WORK_DIR=/scratch/course/2023w300106/juny/day03_GenomeAnnotation/02_augustus augustus run
cd ${WORK_DIR}

### Run samtools faidx to extract the chromosome/scaffold/sequence of interest ###

samtools faidx ${GENOME_DIR}/${GENOME_NAME} ${CHROM_NAME} > ${CHROM_NAME}.fasta

### Augustus runs very slow so it is helpful to divide big reference into chunks, e.g. per chromosome ###
echo "Starting: Augustus . . . (``' 3 ')y "
start_run=$(date +%M)

augustus --species=nematostella_vectensis \
        --strand=both \
        --genemodel=complete \
        --gff3=on \
        --outfile=${WORK_DIR}/${CHROM_NAME}.augustus.gtf \
        ${WORK_DIR}/${CHROM_NAME}.fasta

end_run=$(date +%M)
echo "Augustus finished!"
echo "Elapsed time: $(($end_run-$start_run)) m"

conda deactivate
ml purge