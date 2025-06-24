#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name coverm
#SBATCH --output /scratch/syersin2/Pastobiome_scratch/coverm/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Pastobiome_scratch/coverm/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 50G
#SBATCH --time 02:00:00
#SBATCH --array=1-346

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/CoverM/coverm

# Variables
genome_dir=/scratch/syersin2/Pastobiome_scratch/dereplication_species/dereplicated_genomes
reads_dir=/scratch/syersin2/Pastobiome_scratch/data/bw_cleaned_reads
outdir=/scratch/syersin2/Pastobiome_scratch/coverm/coverm_out
TMPDIR=/scratch/syersin2/Pastobiome_scratch/coverm/tmp

## Array variables
cd ${reads_dir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p ${outdir}/${sample_name}

coverm genome --coupled ${reads_dir}/${sample_name}/${sample_name}_bowtie2_R1_001.fastq.gz \
    ${reads_dir}/${sample_name}/${sample_name}_bowtie2_R2_001.fastq.gz \
    --genome-fasta-directory ${genome_dir} \
    -x fa \
    -m trimmed_mean \
    -p bwa-mem \
    -t 16 \
    -o ${outdir}/${sample_name}/${sample_name}_CoverM.tsv