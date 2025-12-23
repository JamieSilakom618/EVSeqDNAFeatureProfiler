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
# Use relative paths to data folder
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(script_dir)

ev_file_path = os.path.join(repo_root, "data", "gene_fpkm.csv")
srr_file_path = os.path.join(repo_root, "data", "srr5658399_count_with_lenght.xlsx")

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
ev_df["fpkm_evseq"] = (ev_df["reads"] * 1e9) / (ev_df["region_length"] * ev_total)

# SRR RNA-seq FPKM
srr_df["fpkm_srr"] = (srr_df["count"] * 1e9) / (srr_df["region_length"] * srr_total)


# -----------------------------
# 4. Merge tables by gene ID
# -----------------------------
merged_df = pd.merge(
    ev_df[["name", "region_length", "reads", "cpm_evseq", "fpkm_evseq"]],
    srr_df[["gene_id", "count", "cpm_srr", "fpkm_srr"]],
    left_on="name",
    right_on="gene_id",
    how="inner"   # use "outer" to keep all genes
)
#merged_df[merged_df[["cpm_evseq", "cpm_srr"]].isna().any(axis=1)]

# -----------------------------
# 4. Save output
# -----------------------------
output_file = "merged_CPM_table.csv"
data_output = os.path.join(repo_root, "data", output_file)

# Save in both current directory and data folder for flexibility
merged_df.to_csv(output_file, index=False)
merged_df.to_csv(data_output, index=False)

print(f"Merged CPM table saved as {output_file}")
print(f"Also saved to {data_output}")
