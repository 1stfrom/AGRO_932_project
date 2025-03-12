#!/bin/sh
#SBATCH --ntasks-per-node=12
#SBATCH --nodes=1
#SBATCH --mem=68gb
#SBATCH --time=100:00:00
#SBATCH --job-name=map
#SBATCH --array=1-9
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../../log/stdout-%A_%a.log
#SBATCH --output=../../log/stdout-%A_%a.log

cd ../../data/largedata/cycle_30

module load bwa samtools bcftools picard/2.18

SAMPLE_ID=$SLURM_ARRAY_TASK_ID

bwa mem ../ref/Zm-B73-REFERENCE-NAM-5.0.fa T${SAMPLE_ID}_1.fastq T${SAMPLE_ID}_2.fastq | samtools view -bSh - > T${SAMPLE_ID}.bam

samtools sort T${SAMPLE_ID}.bam -o sorted_T${SAMPLE_ID}.bam

samtools index sorted_T${SAMPLE_ID}.bam

picard -Xmx60g MarkDuplicates I=sorted_T${SAMPLE_ID}.bam O=sorted_T${SAMPLE_ID}.picard.bam M=sorted_T${SAMPLE_ID}.picard.matrics AS=true
