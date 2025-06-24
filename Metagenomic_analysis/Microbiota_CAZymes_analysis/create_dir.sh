#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name create_dir
#SBATCH --output /scratch/syersin2/Pastobiome_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Pastobiome_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 10G
#SBATCH --time 00:30:00

cd /scratch/syersin2/Pastobiome_scratch/data/bw_cleaned_reads

for prefix in $(ls *.fastq.gz | sed -E 's/_bowtie2_R[12]_001[.]fastq.gz//' | uniq)
	do
    mkdir ${prefix}
    cp ${prefix}_bowtie2_R1_001.fastq.gz ${prefix}
    cp ${prefix}_bowtie2_R2_001.fastq.gz ${prefix}
done