#!/usr/bin/env python3

import pandas as pd
import sys

def main(file1, file2, out_file, left_col, right_col):
    # Load files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Check columns
    if left_col not in df1.columns:
        sys.exit(f"ERROR: Column '{left_col}' not found in {file1}")
    if right_col not in df2.columns:
        sys.exit(f"ERROR: Column '{right_col}' not found in {file2}")

    # Merge intersecting genes
    merged = pd.merge(
        df1,
        df2,
        left_on=left_col,
        right_on=right_col,
        how="inner"
    )

    # Write output
    merged.to_csv(out_file, index=False)

    print("✅ Intersection completed")
    print(f"File 1: {file1} ({df1.shape[0]} rows)")
    print(f"File 2: {file2} ({df2.shape[0]} rows)")
    print(f"Intersect rows: {merged.shape[0]}")
    print(f"Output: {out_file}")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(
            "Usage:\n"
            "  python intersect_genes_left_right.py "
            "<file1.csv> <file2.csv> <output.csv> <left_gene_col> <right_gene_col>\n\n"
            "Example:\n"
            "  python intersect_genes_left_right.py "
            "gene_99thpercentile_high_abundance.csv "
            "gene_srr_99thpercetle_high_abundance.csv "
            "intersect_genes.csv gene_id Gene"
        )
        sys.exit(1)

    _, file1, file2, out_file, left_col, right_col = sys.argv
    main(file1, file2, out_file, left_col, right_col)

