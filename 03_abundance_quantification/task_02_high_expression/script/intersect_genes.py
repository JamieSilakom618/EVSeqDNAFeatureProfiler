#%%
#!/usr/bin/env python3

import pandas as pd
import sys
import os
import argparse


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TASK_ROOT = os.path.dirname(SCRIPT_DIR)


def resolve_column(df, explicit_col, candidates, file_path):
    if explicit_col:
        if explicit_col in df.columns:
            return explicit_col
        sys.exit(f"ERROR: Column '{explicit_col}' not found in {file_path}. Available: {list(df.columns)}")

    for col in candidates:
        if col in df.columns:
            return col

    sys.exit(
        f"ERROR: Could not auto-detect gene column in {file_path}. "
        f"Available columns: {list(df.columns)}"
    )


def parse_args(argv=None):
    default_file1 = os.path.join(TASK_ROOT, "output_data", "gene_99thpercentile_high_abundance.csv")
    default_file2 = os.path.join(TASK_ROOT, "output_data", "gene_srr_99thpercentile_high_abundance.csv")
    default_out = os.path.join(TASK_ROOT, "output_data", "intersect_genes.csv")

    parser = argparse.ArgumentParser(
        description="Intersect two gene lists by configurable ID columns"
    )
    parser.add_argument("file1", nargs="?", default=default_file1, help="Left CSV file")
    parser.add_argument("file2", nargs="?", default=default_file2, help="Right CSV file")
    parser.add_argument("out_file", nargs="?", default=default_out, help="Output CSV file")
    parser.add_argument("left_col", nargs="?", default=None, help="Left gene ID column")
    parser.add_argument("right_col", nargs="?", default=None, help="Right gene ID column")

    # In notebooks, ignore kernel launcher args (for example: -f <connection.json>)
    if argv is None:
        if "ipykernel" in sys.modules:
            argv = []
        else:
            argv = sys.argv[1:]

    args, unknown = parser.parse_known_args(argv)
    if unknown and "ipykernel" not in sys.modules:
        parser.error(f"unrecognized arguments: {' '.join(unknown)}")

    return args

def main(file1, file2, out_file, left_col, right_col):
    # Load files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Resolve columns (explicit or auto-detected)
    left_col = resolve_column(df1, left_col, ["gene_id", "name", "Gene", "gene", "gene_name"], file1)
    right_col = resolve_column(df2, right_col, ["Gene", "name", "gene_id", "gene", "gene_name"], file2)

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
    print(f"Matched columns: {left_col} (left) vs {right_col} (right)")
    print(f"Intersect rows: {merged.shape[0]}")
    print(f"Output: {out_file}")

if __name__ == "__main__":
    args = parse_args()
    main(args.file1, args.file2, args.out_file, args.left_col, args.right_col)

