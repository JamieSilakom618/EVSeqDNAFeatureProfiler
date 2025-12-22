#!/usr/bin/env python3
"""
FPKM Calculation for EV-seq Analysis

This script calculates Fragments Per Kilobase Million (FPKM) values for genomic features
from read count data for EV-seq analysis.

Author: Nutticha Silakom
Institution: Chulalongkorn University, Bangkok, Thailand
Program: Bioinformatics and Computational Biology, Graduate School
Version: 1.0.0
Date: December 2025

Usage:
    python calculate_fpkm.py

Requirements:
    - Python 3.8+
    - pandas >= 1.3.0
    - samtools (for read count extraction)
    - configparser (standard library)

Input:
    - Read count files in coverage_count/ directory (*.counts format)
    - Aligned BAM file (aligned.mapped.sorted.bam) for read count extraction
    - Configuration parameters from ../config.ini (optional)

Output:
    - FPKM values for each genomic feature class in fpkm_out/ directory
    - Summary statistics for feature classes
"""

import os
import glob
import pandas as pd
import configparser
import sys
import subprocess
from datetime import datetime

# Print script information
print("="*60)
print("FPKM Calculation for EV-seq Analysis")
print("Version 1.0.0")
print(f"Execution time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*60)

# -----------------------------
# Get total mapped reads from BAM file
# -----------------------------
def get_total_mapped_reads(bam_file="aligned.mapped.sorted.bam"):
    """
    Get total mapped reads from BAM file using samtools
    
    Args:
        bam_file (str): Path to BAM file
        
    Returns:
        int: Total number of mapped reads
    """
    try:
        # Run samtools view -c -F 4 to count mapped reads
        result = subprocess.run(
            ["samtools", "view", "-c", "-F", "4", bam_file],
            capture_output=True,
            text=True,
            check=True
        )
        total_reads = int(result.stdout.strip())
        print(f"Total mapped reads from {bam_file}: {total_reads:,}")
        return total_reads
        
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools: {e}")
        print("Make sure samtools is installed and BAM file exists")
        sys.exit(1)
    except ValueError as e:
        print(f"Error parsing samtools output: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: samtools not found. Please install samtools or add it to PATH")
        sys.exit(1)

# Get total mapped reads from BAM file
TOTAL_MAPPED_READS = get_total_mapped_reads()

# Also try to get from config as fallback (for documentation purposes)
config = configparser.ConfigParser()
config_path = '../config.ini'

if os.path.exists(config_path):
    config.read(config_path)
    try:
        config_reads = int(config.get('fpkm_calculation', 'total_mapped_reads', fallback='0'))
        if config_reads > 0:
            print(f"Config file specifies: {config_reads:,} reads")
            if abs(TOTAL_MAPPED_READS - config_reads) > 1000:  # Allow small differences
                print(f"Warning: BAM file read count differs from config by {abs(TOTAL_MAPPED_READS - config_reads):,}")
    except (ValueError, configparser.Error):
        pass

# -----------------------------
# Validate input directories
# -----------------------------
COUNTS_DIR = "coverage_count"
OUT_DIR = "fpkm_out"

if not os.path.exists(COUNTS_DIR):
    print(f"Error: Input directory '{COUNTS_DIR}' not found")
    print("Please ensure read count files are available in the coverage_count directory")
    sys.exit(1)

os.makedirs(OUT_DIR, exist_ok=True)
print(f"Output directory: {OUT_DIR}")

# Check for input files
count_files = glob.glob(os.path.join(COUNTS_DIR, "*.counts"))
if not count_files:
    print(f"Error: No count files (*.counts) found in {COUNTS_DIR}")
    sys.exit(1)
    
print(f"Found {len(count_files)} count files to process")
print(f"Using total mapped reads: {TOTAL_MAPPED_READS:,}")

# -----------------------------
# Store class-level summary
# -----------------------------
class_summary = []

# -----------------------------
# Loop over all *.counts files
# -----------------------------
for counts_file in glob.glob(os.path.join(COUNTS_DIR, "*.counts")):

    feature_class = os.path.basename(counts_file).replace(".counts", "")
    print(f"Processing {feature_class}")

    # Read counts file
    df = pd.read_csv(
        counts_file,
        sep="\t",
        header=None,
        names=[
            "chr",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "reads",
            "bases_covered",
            "region_length",
            "fraction_covered"
        ]
    )

    # Avoid division by zero
    df = df[df["region_length"] > 0].copy()

    # -----------------------------
    # Per-region FPKM
    # -----------------------------
    df["FPKM"] = (
        df["reads"] * 1e9
    ) / (
        df["region_length"] * TOTAL_MAPPED_READS
    )

    # -----------------------------
    # Write per-region output
    # -----------------------------
    out_csv = os.path.join(OUT_DIR, f"{feature_class}_fpkm.csv")
    df_out = df[
        [
            "chr",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "reads",
            "bases_covered",
            "region_length",
            "FPKM"
        ]
    ]

    df_out.to_csv(out_csv, index=False)

    # -----------------------------
    # Feature-class (pooled) FPKM
    # -----------------------------
    total_reads = df["reads"].sum()
    total_length = df["region_length"].sum()

    class_fpkm = (
        total_reads * 1e9
    ) / (
        total_length * TOTAL_MAPPED_READS
    ) if total_length > 0 else 0

    class_summary.append({
        "feature_class": feature_class,
        "total_reads": int(total_reads),
        "total_length": int(total_length),
        "FPKM": class_fpkm
    })

# -----------------------------
# Write class-level summary
# -----------------------------
summary_df = pd.DataFrame(class_summary)
summary_df.to_csv(
    os.path.join(OUT_DIR, "feature_class_FPKM_summary.csv"),
    index=False
)

print("FPKM calculation completed successfully.")

