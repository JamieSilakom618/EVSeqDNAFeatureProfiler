#!/usr/bin/env python3
"""
Prepare Mitochondrial Features for LOLA Region Database

This script processes GFF/GFF3 annotation files to create BED files for mitochondrial 
genomic features, preparing them for LOLA (Locus Overlap Analysis) enrichment testing.

Author: Nutticha Silakom
Institution: Chulalongkorn University, Bangkok, Thailand
Program: Bioinformatics and Computational Biology, Graduate School
Version: 1.0.0
Date: December 2025

Purpose:
Read a GFF/GFF3 file and extract mitochondrial features (chromosome = "Mito") 
to create standardized BED6 files for LOLA region database.

Usage:
    python3 prepare_mtfeatures_forLOLA_regionDB.py input.gff -o output_directory

Input:
    - GFF/GFF3 annotation file (can be gzipped)
    - Output directory path

Output:
    BED6 files for each mitochondrial feature class:
    - gene.bed, tRNA_gene.bed, rRNA_gene.bed, etc.
    Format: chr, start(0-based), end, name, score, strand

Features Processed:
    - Protein-coding genes
    - tRNA genes  
    - rRNA genes
    - Origins of replication
    - Other regulatory elements
"""

from __future__ import annotations
import argparse
import gzip
import os
from collections import defaultdict
from typing import Dict, Iterable, List

# ---- USER MAPPING: keys => list of GFF types to collect ----
# Keys will be used as filenames: <key>.bed
FEATURE_MAPPING: Dict[str, List[str]] = {
    "gene": ["gene"],
    "tRNA_gene": ["tRNA_gene"],
    "rRNA_gene": ["rRNA_gene"],
    "origin_of_replication": [
        "origin_of_replication",
        "ARS",
        "ARS_consensus_sequence"
    ],
    # add more mappings as needed
}
# ----------------------------------------------------------

# optional: additional single-feature names to export as their own file (not needed here)
SINGLE_FEATURES = set()

def open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def parse_attrs(attrstr: str) -> Dict[str, str]:
    """Parse GFF attributes column into a dict. Handles key=value; tolerates '.'"""
    d = {}
    if not attrstr or attrstr.strip() == ".":
        return d
    for part in attrstr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        elif " " in part:
            k, v = part.split(" ", 1)
            d[k] = v.strip('"')
        else:
            d[part] = ""
    return d

def build_reverse_map(mapping: Dict[str, Iterable[str]]) -> Dict[str, List[str]]:
    """Return gff_type -> [mapping_key, ...]"""
    rev = defaultdict(list)
    for key, vals in mapping.items():
        if isinstance(vals, (list, tuple, set)):
            for v in vals:
                rev[str(v)].append(key)
        else:
            rev[str(vals)].append(key)
    return rev

def choose_name(attrd: Dict[str, str], ftype: str, counter: Dict[str, int]) -> str:
    """Choose a name for the BED 'name' column using common attributes"""
    name = attrd.get("Name") or attrd.get("ID") or attrd.get("IDREF") or attrd.get("gene") or ftype
    if not name:
        counter[ftype] += 1
        name = f"{ftype}_{counter[ftype]}"
    return name

def write_beds(outdir: str, beds: Dict[str, List[List[str]]], overwrite: bool = False):
    """Write bed lists (rows are lists of strings) to outdir as <key>.bed"""
    os.makedirs(outdir, exist_ok=True)
    for key, rows in beds.items():
        outpath = os.path.join(outdir, f"{key}.bed")
        mode = "w" if (overwrite or not os.path.exists(outpath)) else "a"
        if mode == "a":
            # if appending and file exists, ensure no duplicate header (we don't write headers anyway)
            pass
        with open(outpath, mode) as outfh:
            for r in rows:
                outfh.write("\t".join(r) + "\n")
        print(f"Wrote {len(rows)} regions to {outpath}")

def extract_mito_features(gff_path: str, mapping: Dict[str, Iterable[str]],
                          single_features: Iterable[str], outdir: str,
                          overwrite: bool = False):
    revmap = build_reverse_map(mapping)
    beds: Dict[str, List[List[str]]] = defaultdict(list)
    counters = defaultdict(int)
    total_in = 0
    total_emitted = 0

    with open_maybe_gz(gff_path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            total_in += 1
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, src, ftype, start_s, end_s, score, strand, phase, attrs = cols[:9]

            # **KEEP ONLY mitochondrial rows where first column EXACTLY equals "Mito"**
            if seqid != "Mito":
                continue

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue

            bed_start = max(0, start - 1)
            bed_end = end

            attrd = parse_attrs(attrs)
            name = choose_name(attrd, ftype, counters)
            score_out = score if score != "." else "0"
            strand_out = strand if strand in {"+", "-", "."} else "."

            # find mapping categories
            cat_keys = revmap.get(ftype, []).copy()
            if ftype in single_features and ftype not in revmap:
                cat_keys.append(ftype)
            if not cat_keys:
                # not mapped -> skip
                continue

            for key in cat_keys:
                # build bed6 row (as strings)
                row = [seqid, str(bed_start), str(bed_end), name, str(score_out), strand_out]
                beds[key].append(row)
                total_emitted += 1

    write_beds(outdir, beds, overwrite=overwrite)
    print(f"Processed {total_in} non-comment lines from {gff_path}")
    print(f"Emitted {total_emitted} mitochondrial regions into {len(beds)} files under {outdir}")

def main():
    p = argparse.ArgumentParser(description="Extract mitochondrial (seqid == 'Mito') features into per-mapping BED files.")
    p.add_argument("gff", help="Input GFF/GFF3 (can be .gz)")
    p.add_argument("-o", "--outdir", default="regionDB", help="Output directory for BED files (default: regionDB)")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing BED files instead of appending")
    args = p.parse_args()

    extract_mito_features(args.gff, FEATURE_MAPPING, SINGLE_FEATURES, args.outdir, overwrite=args.overwrite)

if __name__ == "__main__":
    main()

