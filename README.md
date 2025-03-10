# Literature Source

**Beissinger TM, Hirsch CN, Vaillancourt B, et al.**  
*A genome-wide scan for evidence of selection in a maize population under long-term artificial selection for ear number.*  
Genetics. 2014 Mar;196(3):829-840.  

## Experimental Design

This study conducted a genome-wide selection scan in the *Golden Glow* maize population, which underwent artificial selection for increased ear number over 30 generations. A total of 48 plants were randomly chosen from cycle 0 and cycle 30. The effective population size (\(N_e\)) ranged between **384 and 667 individuals**.

## Materials and Methods

### Germplasm

- The *Golden Glow* maize population has been under selection for increased ear number since **1971**.
- The initial selection intensity was **2.5–5%** in the first 12 cycles and was later increased to **0.5–1%** from the 13th cycle onwards.

### Sequencing Method and Genome Size

- Whole-genome sequencing was performed on 48 randomly selected plants from **cycle 0** and **cycle 30** using **Illumina sequencing technology**.
- Sequencing output:
  - **Cycle 0:** 555,078,520 read pairs from eight sequencing lanes.
  - **Cycle 30:** 652,901,808 read pairs from nine sequencing lanes.

### Ancestral Species

- **Cycle 0** represents the ancestral population.
- **Cycle 30** represents the domesticated population after 30 generations of artificial selection.
- The **reference genome** used for this study was *B73 Version 2*.

## Download Data

To download sequencing data and reference genome files, use the following script:

### **SLURM Batch Script for Data Download**
```bash
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

# Define directories
DATA_DIR30=../data/cycle_30/sra
FASTQ_DIR30=../data/cycle_30/fastq
DATA_DIR0=../data/cycle_0/sra
FASTQ_DIR0=../data/cycle_0/fastq

# SRR IDs for Cycle 30
SRR_IDS30=("SRR801187" "SRR801188" "SRR801189" "SRR801190"
         "SRR801191" "SRR801192" "SRR801194" "SRR801195" "SRR801196")

# Download and convert Cycle 30 data
prefetch --max-size 200G --output-directory $DATA_DIR30 ${SRR_IDS30[@]}

for SRR_ID in ${SRR_IDS30[@]}; do
    fasterq-dump --split-files --threads 16 $DATA_DIR30/$SRR_ID -O $FASTQ_DIR30 &
done

# SRR IDs for Cycle 0
SRR_IDS0=("SRR801164" "SRR801169" "SRR801170" "SRR801174"
         "SRR801175" "SRR801176" "SRR801177" "SRR801178")

# Download and convert Cycle 0 data
prefetch --max-size 200G --output-directory $DATA_DIR0 ${SRR_IDS0[@]}

for SRR_ID in ${SRR_IDS0[@]}; do
    fasterq-dump --split-files --threads 16 $DATA_DIR0/$SRR_ID -O $FASTQ_DIR0 &
done

# Download reference genome
cd ../data/ref
wget -c https://download.maizegdb.org/B73_RefGen_v2/B73_RefGen_v2.fa.gz
gunzip B73_RefGen_v2.fa.gz

# Download reference genome annotation
wget -c https://download.maizegdb.org/B73_RefGen_v2/ZmB73_5a.59_WGS.gff3.gz
gunzip ZmB73_5a.59_WGS.gff3.gz
