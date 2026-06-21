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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORKFLOW_ROOT="$(cd "${TASK_ROOT}/../.." && pwd)"

REGION_DIR="${TASK_ROOT}/input_data/region_for_fpkm"
BAM="${WORKFLOW_ROOT}/aligned.mapped.sorted.bam"
OUTDIR="${TASK_ROOT}/output_data/coverage_count"

# Fallback: also check current working directory
if [[ ! -f "${BAM}" ]]; then
    BAM="aligned.mapped.sorted.bam"
fi

# NOTE: The BAM file (aligned.mapped.sorted.bam) is NOT included in this repository
# because it is too large for version control (typically >10 GB).
# You must supply it before running this script.
if [[ ! -f "${BAM}" ]]; then
    echo "" >&2
    echo "ERROR: BAM file not found: ${BAM}" >&2
    echo "" >&2
    echo "The aligned BAM file is too large to include in this repository." >&2
    echo "Please generate or obtain it and place it at:" >&2
    echo "  ${WORKFLOW_ROOT}/aligned.mapped.sorted.bam" >&2
    echo "" >&2
    exit 1
fi

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

