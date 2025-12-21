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

import sys
import argparse
import requests
import pandas as pd
import time


def parse_args():
    parser = argparse.ArgumentParser(
        description="Lookup yeast genes in SGD using the official locus API"
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
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

    # Load input CSV
    try:
        df = pd.read_csv(args.input)
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

