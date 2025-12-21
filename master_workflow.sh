#!/bin/bash
# master_workflow.sh
# EV-seq Analysis Master Workflow
#
# Author: Nutticha Silakom
# Institution: Chulalongkorn University, Bangkok, Thailand
# Program: Bioinformatics and Computational Biology, Graduate School
# Version: 1.0.0
# 
# This script orchestrates the complete EV-seq analysis pipeline
# Usage: ./master_workflow.sh input_R1.fastq input_R2.fastq sample_name

set -e  # Exit on any error

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_R1.fastq> <input_R2.fastq> <sample_name>"
    echo "Example: $0 SAC11DNA_L2_1.fastq SAC11DNA_L2_2.fastq SAC11DNA"
    exit 1
fi

INPUT_R1=$1
INPUT_R2=$2
SAMPLE_NAME=$3

echo "Starting EV-seq analysis workflow for sample: $SAMPLE_NAME"
echo "=============================================="

# Step 1: Read Processing and Genomic Composition Analysis
echo "Step 1: Read Processing and Genomic Composition Analysis"
cd 01_read_processing_and_genomic_composition
./qc_alignment.sh "$INPUT_R1" "$INPUT_R2" "$SAMPLE_NAME"
cd ..

# Step 2: Feature-Level Locus Enrichment Analysis
echo "Step 2: Feature-Level Locus Enrichment Analysis"
cd 02_locus_enrichment
python prepare_regionBD_forlola.py
python prepare_mtfeatures_forLOLA_regionDB.py
Rscript lola_run.R
cd ..

# Step 3: Quantification of EV DNA Abundance
echo "Step 3: Quantification of EV DNA Abundance"
cd 03_abundance_quantification
./run_bedtools_coverage.sh
python calculate_fpkm.py
python merge_file.py
python correlation_test.py
cd ..

# Collect results
echo "Collecting final results..."
mkdir -p results/$SAMPLE_NAME
cp 01_read_processing_and_genomic_composition/*.html results/$SAMPLE_NAME/
cp 02_locus_enrichment/lola_output/* results/$SAMPLE_NAME/
cp 03_abundance_quantification/fpkm_out/* results/$SAMPLE_NAME/
cp 03_abundance_quantification/merged_CPM_table.csv results/$SAMPLE_NAME/

echo "EV-seq analysis complete!"
echo "Results available in: results/$SAMPLE_NAME/"