#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH -J 10X_gen

# Error report
set -euo pipefail

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Load required tools
echo -e "\n`date` Load modules"
module load bioinfo-tools
module load samtools/1.14
module load sratools/2.10.9
module load cellranger/6.1.2

# Download
# echo -e "\n`date` Start download"
# sam-dump SRR13782529.sra | samtools view -bS - > SRR13782529.bam

# Convert to fastq
echo -e "\n`date` Start conversion to fastq"
mkdir -p $SNIC_TMP/SRR13782529
cd $SNIC_TMP
cellranger bamtofastq --nthreads=10 ${wd}/SRR13782529.bam $SNIC_TMP/SRR13782529

cp -r $SNIC_TMP/SRR13782529 ${wd}

echo -e "\n`date` Done with download"