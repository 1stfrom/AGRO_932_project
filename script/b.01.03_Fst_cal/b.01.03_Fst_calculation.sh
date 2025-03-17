#!/bin/sh
#SBATCH --ntasks-per-node=12
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --time=3:00:00
#SBATCH --job-name=Fst_cal
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../../.err
#SBATCH --output=../../.out

cd ../../cache/SNP/

module load R/4.3
R CMD BATCH ../../script/b.01.03_Fst_cal/b.01.03_Fst_calculation.R
