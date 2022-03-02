#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 16
#SBATCH --mem 64G
#SBATCH -t 72:00:00
#SBATCH -J alignment_counts
#SBATCH --array=1-1%1               # eg 1-56%10 (job 1 to 56, with 10 at a time)

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have

# Error report
set -euo pipefail

# Working directory
wd="$1"

# Data directory
dd="$2"

# Set resources
ref_STAR="$3"

# Load required tools (Detection)
module load TrimGalore/0.6.1
module load star/2.7.9a
module load htseq/0.12.4
module load samtools/1.14

# Assign names for arrays
cd "${wd}"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"

# Generate temporary folder structure
mkdir $SNIC_TMP/processed

# Create processed dir with underlying dirs for sample
cd "$SNIC_TMP/processed/"
mkdir -p "${selected_sample}"

# Run TrimGalore on the fastq
echo -e "\n`date` Filtering and trimming ${selected_sample} ..."
trim_galore --cores 8 \
	--trim-n \
	--fastqc --fastqc_args "--outdir $SNIC_TMP/processed/${selected_sample}/trimgalore/ -t 8" \
	--gzip \
	--paired \
	"${dd}/${selected_sample}_1.fastq.gz" \
	"${dd}/${selected_sample}_2.fastq.gz" \
	-o "$SNIC_TMP/processed/${selected_sample}/trimgalore"

cp "$SNIC_TMP/processed/${selected_sample}/trimgalore/*" "${wd}"

#Run STAR
# echo -e "\n`date` Mapping ${selected_sample} with STAR"
# STAR --genomeDir ${ref_STAR} \
#     --runThreadN 16 \
#     --twopassMode Basic \
#     --readFilesIn \
#     "$SNIC_TMP/processed/${selected_sample}/trimgalore/${selected_sample}_1_val_1.fq.gz" \
#     "$SNIC_TMP/processed/${selected_sample}/trimgalore/${selected_sample}_2_val_2.fq.gz" \
#     --readFilesCommand zcat \
#     --outFileNamePrefix "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}" \
#     --outSAMtype BAM SortedByCoordinate \
#     --outSAMattributes All \
#     --outReadsUnmapped Fastx \
#     --limitBAMsortRAM 150000000000

# # Index by samtools
# echo -e "\n`date` Indexing ${selected_sample} with samtools"
# samtools index -@ 16 "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}Aligned.sortedByCoord.out.bam"

echo -e "\n`date` Finished!"