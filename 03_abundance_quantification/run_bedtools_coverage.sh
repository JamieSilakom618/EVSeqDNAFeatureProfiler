#!/usr/bin/env bash
#
# BEDtools Coverage Analysis for EV-seq Quantification
#
# Author: Nutticha Silakom
# Institution: Chulalongkorn University, Bangkok, Thailand
# Program: Bioinformatics and Computational Biology, Graduate School
# Version: 1.0.0
# Date: December 2025
#
# Purpose:
#   Calculate read coverage for genomic features using BEDtools for FPKM quantification
#   in EV-seq abundance analysis.
#
# Description:
#   - Processes all BED files in region_for_fpkm directory
#   - Calculates coverage of aligned reads for each genomic feature
#   - Generates count files for downstream FPKM calculation
#
# Usage:
#   ./run_bedtools_coverage.sh
#
# Requirements:
#   - BEDtools >= 2.29.0
#   - Aligned BAM file (aligned.mapped.sorted.bam)
#   - BED files defining genomic regions
#
# Output:
#   - Coverage count files in coverage_count/ directory
#   - Format: BED with additional coverage columns
#

set -euo pipefail

REGION_DIR="FPKM/region_for_fpkm"
BAM="aligned.mapped.sorted.bam"
OUTDIR="coverage_count"

mkdir -p "${OUTDIR}"

for bed in "${REGION_DIR}"/*.bed; do
    base=$(basename "${bed}" .bed)
    out="${OUTDIR}/${base}.counts"

    echo "Processing ${bed} → ${out}"

    bedtools coverage \
        -a "${bed}" \
        -b "${BAM}" \
        > "${out}"
done

echo "All coverage files generated successfully."

