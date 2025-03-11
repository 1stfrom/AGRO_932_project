#!/bin/sh
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=60:00:00
#SBATCH --job-name=map
#SBATCH --array=1-12
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../../log/stdout-%A_%a.log
#SBATCH --output=../../log/stdout-%A_%a.log

cd ../../data/largedata/cycle_0

module load bwa samtools bcftools

SAMPLE_ID=$SLURM_ARRAY_TASK_ID

bwa mem ../ref/Zm-B73-REFERENCE-NAM-5.0.fa Z${SAMPLE_ID}_1.fastq Z${SAMPLE_ID}_2.fastq | samtools view -bSh - > Z${SAMPLE_ID}.bam

samtools sort Z${SAMPLE_ID}.bam -o sorted_Z${SAMPLE_ID}.bam

samtools index sorted_Z${SAMPLE_ID}.bam
