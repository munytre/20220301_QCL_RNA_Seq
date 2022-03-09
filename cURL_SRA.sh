#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J download_sra

# Error report
set -euo pipefail

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Download
echo -e "\n`date` Start download"
curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR137/SRR13782529/SRR13782529.sra -o SRR13782529_GSM5106143_02_E4_15_H09_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR137/SRR13782530/SRR13782530.sra -o SRR13782530_GSM5106143_02_E4_15_H09_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR137/SRR13782531/SRR13782531.sra -o SRR13782531_GSM5106143_02_E4_15_H09_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR137/SRR13782532/SRR13782532.sra -o SRR13782532_GSM5106144_03_E3_10_H10_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR137/SRR13782533/SRR13782533.sra -o SRR13782533_GSM5106144_03_E3_10_H10_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR137/SRR13782534/SRR13782534.sra -o SRR13782534_GSM5106144_03_E3_10_H10_Mus_musculus_RNA-Seq.sra

echo -e "\n`date` Done with download"