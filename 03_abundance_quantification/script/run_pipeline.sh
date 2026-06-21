#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python3}"

# ── BAM file pre-flight check ──────────────────────────────────────────────────
# The aligned BAM file is NOT bundled in this repository — it is too large for
# version control (typically >10 GB).  Steps [1/6] and [2/6] require it.
# Place the file at: <project_root>/aligned.mapped.sorted.bam
# See task_01_fpkm/script/run_bedtools_coverage.sh for alignment instructions.
# ────────────────────────────────────────────────────────────────────────────────
BAM="${WORKFLOW_ROOT}/../aligned.mapped.sorted.bam"
if [[ ! -f "${BAM}" ]]; then
    echo "" >&2
    echo "WARNING: aligned.mapped.sorted.bam not found." >&2
    echo "Steps [1/6] and [2/6] will be skipped." >&2
    echo "Run task_01_fpkm/script/run_bedtools_coverage.sh after placing the BAM file." >&2
    echo "" >&2
    SKIP_BAM_STEPS=true
else
    SKIP_BAM_STEPS=false
fi

echo "[1/6] BEDtools coverage"
if [[ "${SKIP_BAM_STEPS}" == true ]]; then
    echo "  → Skipped (BAM file not found)"
else
    bash "${WORKFLOW_ROOT}/task_01_fpkm/script/run_bedtools_coverage.sh"
fi

echo "[2/6] FPKM calculation"
if [[ "${SKIP_BAM_STEPS}" == true ]]; then
    echo "  → Skipped (BAM file not found)"
else
    "${PYTHON_BIN}" "${WORKFLOW_ROOT}/task_01_fpkm/script/calculate_fpkm.py"
fi

echo "[3/6] Merge EV-seq and SRR data"
"${PYTHON_BIN}" "${WORKFLOW_ROOT}/task_02_high_expression/script/merge_file.py"

echo "[4/6] Correlation analysis"
"${PYTHON_BIN}" "${WORKFLOW_ROOT}/task_02_high_expression/script/correlation_test.py"

echo "[5/6] Intersect high-expression gene lists"
"${PYTHON_BIN}" "${WORKFLOW_ROOT}/task_02_high_expression/script/intersect_genes.py" \
  "${WORKFLOW_ROOT}/task_02_high_expression/output_data/gene_99thpercentile_high_abundance.csv" \
  "${WORKFLOW_ROOT}/task_02_high_expression/output_data/gene_srr_99thpercentile_high_abundance.csv" \
  "${WORKFLOW_ROOT}/task_02_high_expression/output_data/intersect_genes.csv" \
  gene_id \
  Gene

echo "[6/6] SGD annotation"
"${PYTHON_BIN}" "${WORKFLOW_ROOT}/task_03_gene_annotation/script/sgd_lookup.py"

echo "Pipeline completed successfully."
