#!/bin/sh
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=36:00:00
#SBATCH --job-name=down30
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../log/30.err
#SBATCH --output=../log/30.out

module load SRAtoolkit/2.11

DATA_DIR=../data/cycle_30/sra
FASTQ_DIR=../data/cycle_30/fastq

SRR_IDS=("SRR801187" "SRR801188" "SRR801189" "SRR801190"
         "SRR801191" "SRR801192" "SRR801194" "SRR801195" "SRR801196")

prefetch --max-size 200G --output-directory $DATA_DIR ${SRR_IDS[@]}

for SRR_ID in ${SRR_IDS[@]}; do
    fasterq-dump --split-files --threads 16 $DATA_DIR/$SRR_ID -O $FASTQ_DIR &
done
wait
echo "All jobs completed!"
