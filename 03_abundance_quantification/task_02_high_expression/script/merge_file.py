#%%
#!/usr/bin/env python3
"""
Merge EV-seq and RNA-seq Data for Correlation Analysis

This script merges EV-seq and RNA-seq count data and calculates CPM (Counts Per Million)
normalization for correlation analysis.

Author: Nutticha Silakom
Institution: Chulalongkorn University, Bangkok, Thailand
Program: Bioinformatics and Computational Biology, Graduate School
Version: 1.0.0
Date: December 2025

Purpose:
1. Load EV-seq and RNA-seq count tables
2. Calculate CPM normalization for both datasets
3. Merge data by gene identifier
4. Prepare merged table for correlation analysis

Usage:
    python merge_file.py

Input:
    - EV-seq FPKM table (data/gene_fpkm.csv)
    - RNA-seq count table (data/srr5658399_count_with_length.xlsx)

Output:
    - merged_CPM_table.csv: Combined CPM and FPKM data for correlation analysis

CPM Calculation:
    CPM = (feature_count / total_mapped_reads) × 1,000,000

FPKM Calculation:
    FPKM = (feature_count × 1,000,000,000) / (gene_length × total_mapped_reads)

Requirements:
    - pandas >= 1.3.0
    - openpyxl (for Excel file reading)
"""

import pandas as pd
import sys
import os
from datetime import datetime

# Print script information
print("="*60)
print("EV-seq and RNA-seq Data Merger")
print("Version 1.0.0")
print(f"Execution time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*60)

# -----------------------------
# 1. Load input count tables
# -----------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
task_root = os.path.dirname(script_dir)
workflow_root = os.path.dirname(task_root)
project_root = os.path.dirname(workflow_root)

ev_candidates = [
    os.path.join(workflow_root, "task_01_fpkm", "output_data", "fpkm_out", "gene_fpkm.csv"),
    os.path.join(project_root, "data", "gene_fpkm.csv"),
]

srr_candidates = [
    os.path.join(task_root, "input_data", "srr5658399_count_with_length.xlsx"),
    os.path.join(task_root, "input_data", "srr5658399_count_with_lenght.xlsx"),
    os.path.join(project_root, "data", "srr5658399_count_with_length.xlsx"),
    os.path.join(project_root, "data", "srr5658399_count_with_lenght.xlsx"),
]

ev_file_path = next((path for path in ev_candidates if os.path.exists(path)), ev_candidates[-1])
srr_file_path = next((path for path in srr_candidates if os.path.exists(path)), srr_candidates[-1])

# Check if files exist
if not os.path.exists(ev_file_path):
    print(f"Error: EV-seq data file not found at {ev_file_path}")
    print("Please place your gene_fpkm.csv file in the data/ folder")
    sys.exit(1)

if not os.path.exists(srr_file_path):
    print(f"Error: RNA-seq data file not found at {srr_file_path}")
    print("Please place your srr5658399_count_with_length.xlsx file in the data/ folder")
    sys.exit(1)

print(f"Loading EV-seq data from: {ev_file_path}")
print(f"Loading RNA-seq data from: {srr_file_path}")

ev_df = pd.read_csv(ev_file_path)     # columns: gene, EV_count
srr_df = pd.read_excel(srr_file_path)      # columns: gene, SRR_exp

# -----------------------------
# 2. Calculate CPM for each dataset
# -----------------------------
# EV-seq CPM
ev_total = ev_df["reads"].sum()
ev_df["cpm_evseq"] = (ev_df["reads"] / ev_total) * 1e6

# SRR RNA-seq CPM
srr_total = srr_df["count"].sum()
srr_df["cpm_srr"] = (srr_df["count"] / srr_total) * 1e6

# -----------------------------
# 3. Calculate FPKM for each dataset
# -----------------------------
# EV-seq FPKM: FPKM = (reads * 1,000,000,000) / (gene_length * total_reads)
ev_df["fpkm_evseq"] = ev_df['FPKM']  
# SRR RNA-seq FPKM
srr_df["fpkm_srr"] = (srr_df["count"] * 1e9) / (srr_df["region_length"] * srr_total)

# -----------------------------
# 3.1 Label top 1% by FPKM (99th percentile cutoff)
# -----------------------------
ev_fpkm_cutoff = ev_df["fpkm_evseq"].quantile(0.99)
srr_fpkm_cutoff = srr_df["fpkm_srr"].quantile(0.99)

ev_df["ev_abundance"] = "not_preferentially_represent"
ev_df.loc[ev_df["fpkm_evseq"] >= ev_fpkm_cutoff, "ev_abundance"] = "preferentially_represented"

srr_df["srr_abundance"] = "not_preferentially_represent"
srr_df.loc[srr_df["fpkm_srr"] >= srr_fpkm_cutoff, "srr_abundance"] = "preferentially_represented"

# -----------------------------
# 4. Merge tables by gene ID
# -----------------------------
merged_df = pd.merge(
    ev_df[["name", "region_length", "reads", "cpm_evseq", "fpkm_evseq", "ev_abundance"]],
    srr_df[["gene_id", "count", "cpm_srr", "fpkm_srr", "srr_abundance"]],
    left_on="name",
    right_on="gene_id",
    how="inner"   # use "outer" to keep all genes
)
#merged_df[merged_df[["cpm_evseq", "cpm_srr"]].isna().any(axis=1)]

# -----------------------------
# 4. Save output
# -----------------------------
output_dir = os.path.join(task_root, "output_data")
os.makedirs(output_dir, exist_ok=True)

output_file = os.path.join(output_dir, "merged_CPM_table.csv")
data_output = os.path.join(project_root, "data", "merged_CPM_table.csv")
ev_top1pct_output = os.path.join(output_dir, "gene_99thpercentile_high_abundance.csv")
srr_top1pct_output = os.path.join(output_dir, "gene_srr_99thpercentile_high_abundance.csv")

# Save in task output folder and legacy data folder for compatibility
merged_df.to_csv(output_file, index=False)
merged_df.to_csv(data_output, index=False)

# Build top-1% tables from merged data for EV-only and SRR-only high expression
ev_top1pct_df = merged_df.loc[
    merged_df["ev_abundance"] == "preferentially_represented",
    ["name", "region_length", "reads", "cpm_evseq", "fpkm_evseq", "ev_abundance"],
].copy()
ev_top1pct_df.insert(0, "gene_id", ev_top1pct_df["name"])

srr_top1pct_df = merged_df.loc[
    merged_df["srr_abundance"] == "preferentially_represented",
    ["name", "region_length", "reads", "count", "cpm_srr", "fpkm_srr", "srr_abundance"],
].copy()
srr_top1pct_df.insert(0, "Gene", srr_top1pct_df["name"])

ev_top1pct_df.to_csv(ev_top1pct_output, index=False)
srr_top1pct_df.to_csv(srr_top1pct_output, index=False)

print(f"Merged CPM table saved as {output_file}")
print(f"Also saved to {data_output}")
print(f"EV top 1% table saved as {ev_top1pct_output} ({len(ev_top1pct_df)} rows)")
print(f"SRR top 1% table saved as {srr_top1pct_output} ({len(srr_top1pct_df)} rows)")
