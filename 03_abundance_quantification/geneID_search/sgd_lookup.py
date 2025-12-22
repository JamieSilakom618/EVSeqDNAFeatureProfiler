#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SGD Gene Lookup Tool for EV-seq Analysis

This script queries the Saccharomyces Genome Database (SGD) to retrieve gene annotations
for gene identifiers in EV-seq analysis results.

Author: Nutticha Silakom
Institution: Chulalongkorn University, Bangkok, Thailand
Program: Bioinformatics and Computational Biology, Graduate School
Version: 1.0.0
Date: December 2025

Purpose:
Query SGD (Saccharomyces Genome Database) API to annotate gene lists with:
- Standard gene names
- Gene descriptions
- Systematic names
- Functional annotations

Usage:
    python sgd_lookup.py -i input.csv -o output.csv -c gene_column

Examples:
    # Basic usage with default files
    python sgd_lookup.py

    # Custom input and output files
    python sgd_lookup.py -i ../data/genes.csv -o genes_with_sgd_annotation.csv -c gene_name
    
    # Using gene IDs from different column
    python sgd_lookup.py -i gene_list.csv -o annotated_genes.csv -c systematic_name

Input:
    - CSV file with gene identifiers
    - Column name containing gene IDs

Output:
    - Annotated CSV with SGD gene information

API:
    Uses SGD REST API (https://www.yeastgenome.org/backend/locus/)

Requirements:
    - requests >= 2.25.0
    - pandas >= 1.3.0
    - Internet connection for SGD API access
"""
import os
import sys
import argparse
import requests
import pandas as pd
import time


def parse_args():
    parser = argparse.ArgumentParser(
        description="Lookup yeast genes in SGD using the official locus API"
    )
    parser.add_argument("-i", "--input", 
                       default="../fpkm_out/gene_mt_fpkm.csv",
                       help="Input CSV file (default: ../fpkm_out/gene_mt_fpkm.csv)")
    parser.add_argument("-o", "--output", 
                       default="gene_mt_fpkm_w_name.csv",
                       help="Output CSV file (default: gene_mt_fpkm_w_name.csv)")
    parser.add_argument(
        "-c", "--column",
        default="name",
        help="Column containing gene names (default: name)"
    )
    return parser.parse_args()


def query_sgd_locus(locus):
    """Query SGD locus API for a single gene"""
    url = f"https://www.yeastgenome.org/backend/locus/{locus}"
    headers = {"accept": "application/json"}

    r = requests.get(url, headers=headers, timeout=30)
    if r.status_code != 200:
        return None

    data = r.json()

    # ---- robust description handling ----
    desc = None
    description_field = data.get("description")

    if isinstance(description_field, str):
        desc = description_field
    elif isinstance(description_field, dict):
        refs = description_field.get("references", [])
        if refs and isinstance(refs, list):
            desc = refs[0].get("citation")

    return {
        "query_name": locus,
        "sgdid": data.get("sgdid"),
        "systematic_name": data.get("systematic_name"),
        "gene_name": data.get("gene_name"),
        "format_name": data.get("format_name"),
        "description": desc
    }


def main():
    args = parse_args()

    # Check if input file exists, try data folder as backup
    input_file = args.input
    if not os.path.exists(input_file):
        # Try relative path from data folder
        backup_path = f"../../data/{os.path.basename(input_file)}"
        if os.path.exists(backup_path):
            input_file = backup_path
            print(f"Using backup file from data: {input_file}")
        else:
            print(f"Error: Input file not found: {args.input}")
            print(f"Also checked: {backup_path}")
            sys.exit(1)

    # Load input CSV
    try:
        df = pd.read_csv(input_file)
        print(f"Loaded input file: {input_file}")
    except Exception as e:
        sys.exit(f"ERROR reading input CSV: {e}")

    if args.column not in df.columns:
        sys.exit(f"ERROR: column '{args.column}' not found")

    genes = (
        df[args.column]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    results = []
    not_found = []

    for gene in genes:
        print(f"Querying SGD: {gene}", file=sys.stderr)
        res = query_sgd_locus(gene)
        if res:
            results.append(res)
        else:
            not_found.append(gene)

        # be polite to SGD servers
        time.sleep(0.2)

    df_out = pd.DataFrame(results)
    df_out.to_csv(args.output, index=False)

    print(f"\nSaved: {args.output}")
    print(f"Found: {len(results)} genes")
    print(f"Not found: {len(not_found)}")

    if not_found:
        with open("sgd_not_found.txt", "w") as f:
            f.write("\n".join(not_found))
        print("Missing genes written to: sgd_not_found.txt")


if __name__ == "__main__":
    main()

