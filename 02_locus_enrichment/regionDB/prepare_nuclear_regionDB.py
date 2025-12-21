#!/usr/bin/env python3
"""
Prepare Nuclear Region Database for LOLA Analysis

This script creates standardized BED files from GFF annotation for nuclear chromosomes,
excluding mitochondrial features, to support LOLA enrichment analysis.

Author: Nutticha Silakom
Institution: Chulalongkorn University, Bangkok, Thailand
Program: Bioinformatics and Computational Biology, Graduate School
Version: 1.0.0
Date: December 2025

Purpose:
Create separate BED files for each genomic feature class from GFF annotation,
specifically for nuclear chromosomes (excluding mitochondrial DNA).

Processing:
- Reads GFF/GFF3 annotation file
- Converts coordinates: GFF (1-based inclusive) → BED (0-based half-open)
- Filters out mitochondrial chromosomes
- Creates separate BED files for each feature type

Input:
    - GFF annotation file (s288c_annotation_genome.gff)

Output:
    - Separate BED files for each feature class in chr_regions/
    - Format: BED6 (chr, start, end, name, score, strand)

Feature Classes Processed:
    - Genes, pseudogenes, centromeres
    - Non-coding RNAs (rRNA, tRNA, ncRNA)
    - Regulatory regions and origins of replication
    - Mobile genetic elements and LTRs
"""

import re
import os
from collections import defaultdict

# Hard-coded paths
GFF_PATH = "/home/nutticha/exo_seq/s288c_annotation_genome.gff"
OUTDIR = "./chr_regions/"

FEATURE_MAP = {
    "gene": "gene.bed",
    "pseudogene": "pseudogene.bed",
    "centromere": "centromere.bed",
    "long_terminal_repeat": "long_terminal_repeat.bed",
    "mobile_genetic_element": "mobile_genetic_element.bed",
    "ncRNA": "ncRNA.bed",
    "origin_of_replication": "origin_of_replication.bed",
    "regulatory_region": "regulatory_region.bed",
    "rRNA": "rRNA.bed",
    "tRNA": "tRNA.bed",
    "snRNA": "snRNA.bed",
    "snoRNA": "sno.bed",
    "telomere": "telomere.bed",
}

# regex to detect mitochondrial chromosomes (case-insensitive)
PAT_MITO = re.compile(r"Mito|mito|mt|chrM|mtdna|mitochond", re.I)


def parse_gff(gff_path, feature_keys):
    """
    Parse GFF and return dict: feature -> list of (chrom, start0, end)
    """
    beds = defaultdict(list)

    with open(gff_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            chrom = cols[0].strip()
            feature = cols[2].strip()

            # skip mitochondrial contigs
            if PAT_MITO.search(chrom):
                continue

            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue

            if feature in feature_keys:
                start0 = start - 1
                beds[feature].append((chrom, start0, end))

    return beds


def sort_bed(entries):
    """Sort entries by chromosome and start coordinate."""
    def chrom_key(c):
        m = re.match(r"^(?:chr)?(\d+)$", c)
        if m:
            return (0, int(m.group(1)))
        return (1, c)

    return sorted(entries, key=lambda x: (chrom_key(x[0]), x[1], x[2]))


def write_bed(path, entries):
    """Write sorted BED."""
    with open(path, "w") as out:
        for chrom, start, end in entries:
            out.write(f"{chrom}\t{start}\t{end}\n")


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    features = set(FEATURE_MAP.keys())
    beds = parse_gff(GFF_PATH, features)

    for feat in sorted(features):
        rows = beds.get(feat, [])
        rows_sorted = sort_bed(rows)
        out_name = FEATURE_MAP[feat]
        out_path = os.path.join(OUTDIR, out_name)

        write_bed(out_path, rows_sorted)
        print(f"Saved {out_path}   (feature={feat}, n={len(rows_sorted)})")


if __name__ == "__main__":
    main()

