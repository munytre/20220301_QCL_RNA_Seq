#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 72:00:00
#SBATCH -J alignment_counts

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have in shell

# Error report
set -euo pipefail

# Working directory
wd="$1"

# Set resources
ref_STAR="$2"
ref_gtf="$3"

# Load required tools (Detection)
module load bioinfo-tools
module load TrimGalore/0.6.1
module load star/2.7.9a
module load samtools/1.14
module load htseq/0.12.4
module load GATK/4.2.0.0
module load picard/2.23.4
# java/sun_jdk1.8.0_151
# Python 3.7.2
# CutAdapt 2.3
# FastQC 0.11.8

# Assign names for arrays
cd "${wd}"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"
links=($(cat ftp_links.txt))
ftp_link=${links[$((SLURM_ARRAY_TASK_ID-1))]}

# Generate temporary folder structure
mkdir -p $SNIC_TMP/{processed,data}

# Download data
cd "$SNIC_TMP/data/"
mkdir -p ${selected_sample}
cd ${selected_sample}
wget "${ftp_link}${selected_sample}/${selected_sample}_1.fastq.gz"
wget "${ftp_link}${selected_sample}/${selected_sample}_2.fastq.gz"

# Create processed dir with underlying dirs for sample
cd "$SNIC_TMP/processed/"
mkdir -p "${selected_sample}"

# Run TrimGalore on the fastq
echo -e "\n`date` Filtering and trimming ${selected_sample} ..."
trim_galore --cores 8 \
	--trim-n \
	--fastqc \
	--gzip \
	--paired \
	"$SNIC_TMP/data/${selected_sample}/${selected_sample}_1.fastq.gz" \
	"$SNIC_TMP/data/${selected_sample}/${selected_sample}_2.fastq.gz" \
	-o "$SNIC_TMP/processed/${selected_sample}/trimgalore"

mkdir -p ${wd}/processed/${selected_sample}/trimgalore/
cp $SNIC_TMP/processed/${selected_sample}/trimgalore/*.html ${wd}/processed/${selected_sample}/trimgalore/

#Run STAR
echo -e "\n`date` Mapping ${selected_sample} with STAR"
STAR --genomeDir ${ref_STAR} \
    --runThreadN 16 \
    --twopassMode Basic \
    --readFilesIn \
    "$SNIC_TMP/processed/${selected_sample}/trimgalore/${selected_sample}_1_val_1.fq.gz" \
    "$SNIC_TMP/processed/${selected_sample}/trimgalore/${selected_sample}_2_val_2.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --limitBAMsortRAM 100000000000

# Index by samtools
echo -e "\n`date` Indexing ${selected_sample} with samtools"
samtools index -@ 16 "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}Aligned.sortedByCoord.out.bam"

mkdir -p ${wd}/processed/${selected_sample}/STAR/
cp -R $SNIC_TMP/processed/${selected_sample}/STAR/* ${wd}/processed/${selected_sample}/STAR

# Counting by htseq-count
echo -e "\n`date` Counting raw expression values of ${selected_sample} with htseq-count"
mkdir -p ${wd}/processed/${selected_sample}/htseq_count/
htseq-count -n 16 \
    --order pos \
    --stranded reverse \
    --counts_output "${wd}/processed/${selected_sample}/htseq_count/${selected_sample}.count" \
    "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}Aligned.sortedByCoord.out.bam" \
    "${ref_gtf}"

java -jar $PICARD_ROOT/picard.jar MarkDuplicates

echo -e "\n`date` Finished!"