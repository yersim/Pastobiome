#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name drep_species
#SBATCH --output /scratch/syersin2/Pastobiome_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Pastobiome_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 50G
#SBATCH --time 08:00:00

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/dRep/drep

# Variables
drep_dir=/scratch/syersin2/Pastobiome_scratch/dereplication_species
genomes=/scratch/syersin2/Pastobiome_scratch/data/MAGs
genome_info=/scratch/syersin2/Pastobiome_scratch/data/tables

# Dereplicate 
# Using default for pc = 0.9, sc = 0.95, and nc=0.1
dRep dereplicate ${drep_dir} \
    -g ${genomes}/*.fa \
    --genomeInfo ${genome_info}/YERS23-2.checkm_dRep.csv \
    -p 32 \
    --S_algorithm fastANI \
    -comp 75 -con 10