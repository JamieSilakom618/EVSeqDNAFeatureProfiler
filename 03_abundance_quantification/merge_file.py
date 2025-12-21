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
    - EV-seq FPKM table (gene_fpkm.csv)
    - RNA-seq count table (Excel format)

Output:
    - merged_CPM_table.csv: Combined CPM data for correlation analysis

CPM Calculation:
    CPM = (feature_count / total_mapped_reads) × 1,000,000

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
print("Journal of Extracellular Vesicles - Supporting Code")
print("Version 1.0.0")
print(f"Execution time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*60)

# -----------------------------
# 1. Load input count tables
# -----------------------------
ev_df = pd.read_csv("/Users/nuttichasilakom/Desktop/Paper_projects/natapols_tasks/EVSeq/exo_seq/FPKM/fpkm_out/gene_fpkm.csv")     # columns: gene, EV_count
srr_df = pd.read_excel("/Users/nuttichasilakom/Desktop/Paper_projects/natapols_tasks/EVSeq/exo_seq/FPKM/find lenght/srr5658399_count_with_lenght.xlsx")      # columns: gene, SRR_exp

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
# 3. Merge tables by gene ID
# -----------------------------
merged_df = pd.merge(
    ev_df[["name", "reads", "cpm_evseq"]],
    srr_df[["gene_id", "count", "cpm_srr"]],
    left_on="name",
    right_on="gene_id",
    how="inner"   # use "outer" to keep all genes
)
#merged_df[merged_df[["cpm_evseq", "cpm_srr"]].isna().any(axis=1)]

# -----------------------------
# 4. Save output
# -----------------------------
merged_df.to_csv("merged_CPM_table.csv", index=False)

print("Merged CPM table saved as merged_CPM_table.csv")
