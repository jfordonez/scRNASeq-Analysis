#!/bin/bash
#
#SBATCH --job-name=fannotation
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --partition=basic
#SBATCH --output=log/annot-%j.out
#SBATCH --error=log/annot-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=mail@univie.ac.at
# #SBATCH --test-only


## SUMMARY: This script performs the annotation of the transcriptome from the genome
## (generated based on the gff information)

##load modules

module load gffread #used to extract specific data drom the genome guided by the gff file
module load transdecoder/5.7.1
module load interproscan/5.65-97.0-11.0.4
module load conda/miniconda3
conda activate eggnog-mapper-2.1.12
module list
echo "Conda env activated: $CONDA_DEFAULT_ENV"

#go to the preferred working directory
WORK_DIR=/scratch/course/2023w300106/juny/day05/01_Annotation
cd ${WORK_DIR}

#### Set the paths and variables (if any) ###
GENOME_DIR=/scratch/course/2023w300106/juny/day01/01_dataset/ncbi_dataset/data/GCF_932526225.1 #path to the directory containing the fastq files
GFF_DIR=/scratch/course/2023w300106/juny/day04/01_TPM
OUT_DIR=/scratch/course/2023w300106/juny/day05/01_Annotation

RUN_NAME=Nvec_SRR24348406 ##Input run name
GENOME_FILE=GCF_932526225.1_jaNemVect1.1_genomic.fna
GFF_FILE=genomic.gff

##----------------------------------------------------------------------------------------------------------##
## 1. Generate the transcriptome from the genome based on the gff information ##

echo "Starting: gffread . . . (o,,w,,o ) "

gffread -g ${GENOME_DIR}/${GENOME_FILE} \
      -w ${RUN_NAME}_transcipts.fa \
      ${GFF_DIR}/${GFF_FILE}

echo "Run ended.."
echo ${DATE}

##----------------------------------------------------------------------------------------------------------##
## 2. Run Transdecoder  to predict the amino acid sequences, which will be used in the succeeding analyses ##
start_run=$(date +%M)

echo "Starting: Transdecoder run . . ."

TransDecoder.LongOrfs -S \
      -t ${RUN_NAME}_transcipts.fa \
      --complete_orfs_only

TransDecoder.Predict -t ${RUN_NAME}_transcipts.fa \
      --single_best_only

# remove ".p1" ending in the transcript header ID and the * stop codon to  make it easily readable
# in downstream analyses

sed -e 's/\.p[0-9]\+/ /' -e 's/*//' ${RUN_NAME}_transcripts.fa.transdecoder.pep \
      > ${RUN_NAME}_prot_transdecoder.fa

echo "Transdecoder ended.."

end_run=$(date +%M)
echo "Elapsed time: $(($end_run-$start_run)) m"


##--------------------------------------------------------------------------------------------------------------------##
## 3. Perform Interproscan (IPS) to provide annotation based on search against several protein signature databases
## This will provide more information on what the genes are and their function (based on the data of classical model organisms)
start_run=$(date +%M)

echo "Start interpro scan . . ."

## To have a tidy directory, make new dir to put all IPS result
mkdir ${OUT_DIR}/ips

interproscan.sh -b ${OUT_DIR}/ips \
      -cpu ${SLURM_CPUS_PER_TASK} \
      -etra \
      -goterms \
        # -f TSV,GFF3 \
      -i ${RUN_NAME}_prot_transdecoder.fasta \
      -iprlookup \
      -pa \
      -seqtype p

echo "Interproscan ended.."
end_run=$(date +%M)
echo "Elapsed time: $(($end_run-$start_run)) m"

##----------------------------------------------------------------------------------------------------------##
## 4. Another method of annotating genes is searching against orthologous sequences from eggNOG database

## Tidying
mkdir emapper_out
OUTDIR_EGGNOG=emapper_out

echo "Starting: eggnog mapper. . ."
start_run=$(date +%M)

emapper.py --cpu ${SLURM_CPUS_PER_TASK} \
        -i ${RUN_NAME}_prot_transdecoder.fasta \
        --itype proteins \
        --data_dir /scratch/mirror/eggnog-mapper/2.1.12 \
        --pident 60 \
        --query_cover 60 \
        --subject_cover 70 \
        --tax_scope 'eukaryota_broad' \
        --go_evidence all \
        --output ${RUN_NAME}_emapper-out \
        --output_dir ${OUTDIR_EGGNOG}
        --temp_dir ${TMPDIR} \
        --override

echo "eggnog mapper ended.."

end_run=$(date +%M)
echo "Elapsed time: $(($end_run-$start_run)) s"

conda deactivate
ml purge