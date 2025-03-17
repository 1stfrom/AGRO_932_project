#!/bin/sh
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --job-name=merge
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../../log/.err
#SBATCH --output=../../log/.out

cd ../../cache/
module load bcftools
bcftools merge -Ob -o merged.bcf ../data/largedata/cycle_30/T{1..9}.bcf ../data/largedata/cycle_0/Z{1..8}.bcf
bcftools index merged.bcf
# get a summary of the VCF file
bcftools stats merged.bcf > bcf_summary.txt

# Filter for Biallelic SNPs
bcftools view merged.bcf -m2 -M2 -v snps -Ob -o merged_biallelic_snps.bcf

# get a summary again
bcftools stats  merged_biallelic_snps.bcf > bcf_sum2.txt

# convert to text file
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP[\t%GT]\n' merged_biallelic_snps.bcf > snp_calls.txt
