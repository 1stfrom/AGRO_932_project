#!/bin/sh
#SBATCH --ntasks-per-node=6
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --job-name=indexref
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../log/indexref.err
#SBATCH --output=../log/indexref.out

cd ../data/largedata/ref

module load bwa samtools bcftools


bwa index Zm-B73-REFERENCE-NAM-5.0.fa


