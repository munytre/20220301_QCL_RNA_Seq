#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J 10X_gen

# Error report
set -euo pipefail

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Download
echo -e "\n`date` Start download"
sam-dump SRR13782529 | samtools view -bS - > SRR13782529.bam

echo -e "\n`date` Done with download"