#!/usr/bin/env bash
set -euo pipefail

############################################################
# EV-seq Step 1: Quality Control and Alignment
#
# Author: Nutticha Silakom
# Institution: Chulalongkorn University, Bangkok, Thailand
# Program: Bioinformatics and Computational Biology, Graduate School
# Version: 1.0.0
# Date: December 2025
#
# Description:
#   Performs quality control and genome alignment for EV-seq data:
#   - Run FastQC on paired-end FASTQ files
#   - Align reads to S. cerevisiae reference genome (S288C)
#   - Filter for primary mapped alignments only
#   - Generate sorted and indexed BAM files
#   - Output basic mapping statistics and genomic composition
#
# Usage:
#   ./qc_alignment.sh input_R1.fastq input_R2.fastq output_prefix
#
# Requirements:
#   - FastQC, BWA-MEM, SAMtools
#   - S288C reference genome indexed with BWA
#
# Filtering Strategy:
#   Exclude SAM flags:
#     4    unmapped reads
#     256  secondary alignments  
#     2048 supplementary alignments
#   Total exclusion flag: -F 2308
#
# Output:
#   - FastQC reports (.html)
#   - Sorted BAM file with index
#   - Alignment statistics
#   - Genomic composition metrics
#
# Notes:
#   - Raw FASTQ and BAM files are not committed to GitHub due to the size limited
#
# Usage:
#   bash qc_alignment.sh <sample_id> <fastq_r1> <fastq_r2> <reference_fasta> [threads]
#
# Example:
#   bash qc_alignment.sh \
#     sample01 \
#     data/raw/sample01_R1.fq.gz \
#     data/raw/sample01_R2.fq.gz \
#     reference/genome.fa \
#     8
############################################################

# -----------------------------
# Parse arguments
# -----------------------------
if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <sample_id> <fastq_r1> <fastq_r2> <reference_fasta> [threads]"
  exit 1
fi

SAMPLE_ID="$1"
FASTQ_R1="$2"
FASTQ_R2="$3"
REF_FASTA="$4"
THREADS="${5:-4}"

# -----------------------------
# Output directories
# -----------------------------
QC_DIR="results/fastqc_out"
ALIGN_DIR="results/alignment_out"

mkdir -p "$QC_DIR" "$ALIGN_DIR"

# -----------------------------
# Step 1: FastQC
# -----------------------------
echo "[INFO] Running FastQC for sample: ${SAMPLE_ID}"

fastqc \
  -t "$THREADS" \
  -o "$QC_DIR" \
  "$FASTQ_R1" "$FASTQ_R2"

# -----------------------------
# Step 2: Alignment (BWA-MEM)
# -----------------------------
echo "[INFO] Aligning reads to reference genome"

bwa mem \
  -t "$THREADS" \
  "$REF_FASTA" \
  "$FASTQ_R1" "$FASTQ_R2" \
  > "$ALIGN_DIR/${SAMPLE_ID}.sam"


# -----------------------------
# Step 3: Filter + sort BAM
#   Keep only primary mapped alignments
# -----------------------------
echo "[INFO] Filtering (SAM flag -F 2308), sorting, and indexing BAM"

samtools view -@ "$THREADS" -b -F 2308 \
  "$ALIGN_DIR/${SAMPLE_ID}.sam" \
  | samtools sort -@ "$THREADS" \
    -o "$ALIGN_DIR/${SAMPLE_ID}.sorted.bam"

samtools index "$ALIGN_DIR/${SAMPLE_ID}.sorted.bam"

# -----------------------------
# Step 4: Mapping statistics
# -----------------------------
echo "[INFO] Generating mapping statistics"

samtools view -c -F 4 \
  "$ALIGN_DIR/${SAMPLE_ID}.sorted.bam" \
  > "$ALIGN_DIR/${SAMPLE_ID}.mapped_read_count.txt"

samtools idxstats \
  "$ALIGN_DIR/${SAMPLE_ID}.sorted.bam" \
  > "$ALIGN_DIR/${SAMPLE_ID}.idxstats.txt"

samtools flagstat \
  "$ALIGN_DIR/${SAMPLE_ID}.sorted.bam" \
  > "$ALIGN_DIR/${SAMPLE_ID}.flagstat.txt"

# -----------------------------
# Cleanup
# -----------------------------
rm "$ALIGN_DIR/${SAMPLE_ID}.sam"

echo "[DONE] QC and alignment completed for sample: ${SAMPLE_ID}"

