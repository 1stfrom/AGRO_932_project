#!/bin/sh
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=36:00:00
#SBATCH --job-name=downloadseq
#SBATCH --mail-user=nathanchu@huskers.unl.edu
#SBATCH --mail-type=ALL
#SBATCH --error=../log/.err
#SBATCH --output=../log/.out

module load SRAtoolkit/2.11

DATA_DIR=../data/cycle_0/sra
FASTQ_DIR=../data/cycle_0/fastq

SRR_IDS=("SRR801164" "SRR801169" "SRR801170" "SRR801174"
         "SRR801175" "SRR801176" "SRR801177" "SRR801178")

prefetch --max-size 200G --output-directory $DATA_DIR ${SRR_IDS[@]}

for SRR_ID in ${SRR_IDS[@]}; do
    fasterq-dump --split-files --threads 16 $DATA_DIR/$SRR_ID -O $FASTQ_DIR &
done
wait
echo "All jobs completed!"
