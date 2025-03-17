#!/bin/sh
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --job-name=SNP
#SBATCH --array=1-8
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../../log/stdout-%A_%a.log
#SBATCH --output=../../log/stdout-%A_%a.log

cd ../../data/largedata/cycle_0
module load samtools bcftools
SAMPLE_ID=$SLURM_ARRAY_TASK_ID

# Step 1: Index the sorted BAM file (if not already indexed)
if [ ! -f sorted_Z${SAMPLE_ID}.bam.bai ]; then
    samtools index sorted_Z${SAMPLE_ID}.bam
fi

# Step 2: Generate mpileup file
bcftools mpileup -f ../ref/Zm-B73-REFERENCE-NAM-5.0.fa sorted_Z${SAMPLE_ID}.bam > Z${SAMPLE_ID}.mpileup

# Step 3: Call SNPs with conserved calling and output everything
bcftools call -c -Ob -o Z${SAMPLE_ID}.bcf Z${SAMPLE_ID}.mpileup

# Step 4: Index the VCF file
bcftools index Z${SAMPLE_ID}.bcf
