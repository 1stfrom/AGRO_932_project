#!/bin/sh
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=60:00:00
#SBATCH --job-name=SNP
#SBATCH --array=1-12
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../../log/stdout-%A_%a.log
#SBATCH --output=../../log/stdout-%A_%a.log

cd ../../data/largedata/cycle_30
module load samtools bcftools
SAMPLE_ID=$SLURM_ARRAY_TASK_ID

# Step 1: Index the sorted BAM file (if not already indexed)
if [ ! -f sorted_T${SAMPLE_ID}.bam.bai ]; then
    samtools index sorted_T${SAMPLE_ID}.bam
fi

# Step 2: Generate mpileup file
bcftools mpileup -f ../ref/Zm-B73-REFERENCE-NAM-5.0.fa sorted_T${SAMPLE_ID}.bam > T${SAMPLE_ID}.mpileup

# Step 3: Call SNPs with conserved calling and output everything
bcftools call -c -Ob -o T${SAMPLE_ID}.bcf T${SAMPLE_ID}.mpileup

# Step 4: Index the VCF file
bcftools index T${SAMPLE_ID}.bcf
