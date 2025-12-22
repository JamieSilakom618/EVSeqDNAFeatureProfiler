#!/usr/bin/env python3
"""
Prepare Nuclear Genomic Features for LOLA Region Database

This script converts GFF annotation files to BED format for nuclear genomic features,
creating region databases suitable for LOLA (Locus Overlap Analysis) enrichment testing.

Author: Nutticha Silakom
Institution: Chulalongkorn University, Bangkok, Thailand
Program: Bioinformatics and Computational Biology, Graduate School
Version: 1.0.0
Date: December 2025

Purpose:
Convert GFF to standardized BED files for nuclear chromosomes, excluding mitochondrial
features. Creates separate BED files for each feature class for LOLA analysis.

Usage:
    python3 prepare_regionBD_forlola.py input.gff -o output_directory

Processing:
- Feature mapping: Maps GFF feature types to standardized feature classes
- Coordinate conversion: GFF (1-based inclusive) → BED (0-based, half-open)
- Output format: BED6 (chrom, start, end, name, score, strand)
- Mitochondrial contigs are excluded from nuclear analysis

Input:
    - GFF/GFF3 annotation file
    - Output directory path

Output:
    BED6 files for nuclear feature classes:
    - gene.bed, ncRNA_gene.bed, pseudogene.bed, etc.
    
Feature Classes:
    - Protein-coding genes
    - Non-coding RNAs (rRNA, tRNA, snoRNA, snRNA, ncRNA)
    - Regulatory elements (centromeres, telomeres, origins)
    - Transposable elements
    - Pseudogenes
"""

import sys
import os
import argparse
from collections import defaultdict

# ---- User mapping: change as needed ----
# This reflects the structure you provided: keys are output names, values are lists of GFF types to match.
FEATURE_MAPPING = {"chromosome": ["chromosome"]}
#FEATURE_MAPPING = {
#    "gene": ["gene"],
#    "pseudogene": ["pseudogene"],
#    "ncRNA_gene": ["ncRNA_gene"],
#    "snoRNA_gene": ["snoRNA_gene"],
#    "snRNA_gene": ["snRNA_gene"],
#    "tRNA_gene": ["tRNA_gene"],
#    "rRNA_gene": ["rRNA_gene"],
#    "Replication_origins": [
#        "origin_of_replication",
#        "ARS",
#        "ARS_consensus_sequence"
#    ],
#    "transposable_elements": [
#        "transposable_element_gene",
#        "LTR_retrotransposon",
#        "long_terminal_repeat"
#    ],
#    "Mating_loci": [
#        "mating_type_region",
#        "silent_mating_type_cassette_array"
#    ],
#    "telomere": ["telomere", "telomeric_repeat"],
#    "centromere": ["centromere", "centromere_DNA_Element_I", "centromere_DNA_Element_II", "centromere_DNA_Element_III"],
    # add more mappings if needed
#}

# Also accept single-feature names that should be exported under their own file.
# If you prefer all mapping in FEATURE_MAPPING as above, you can leave this empty.
SINGLE_FEATURES = set([
    # leave empty if you included them in FEATURE_MAPPING already
])

# ---- end mapping ----

def is_mito(seqid: str) -> bool:
    """Heuristic to decide if the seqid is mitochondrial. Adjust if needed for your assembly names."""
    if not seqid:
        return False
    s = seqid.lower()
    # common patterns: 'mito', 'mt', 'chrm', 'chrmt', 'mitochondrion'
    if s.startswith("mito") or s.startswith("mt") or s.startswith("chrm") or "mitochond" in s:
        return True
    # sometimes assemblies use 'chrM' or 'M'
    if s == "m" or s == "chrm" or s == "chrm_m":
        return True
    return False

def parse_attrs(attrstr):
    """Parse GFF attribute column into dict (supports key=value; tolerates '.' )"""
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
            d[k] = v
        else:
            d[part] = ""
    return d

def build_reverse_map(feature_mapping):
    """
    Build dict: feature_type -> list of category_keys
    So if a GFF type matches multiple categories, it will be placed into each category.
    """
    rev = defaultdict(list)
    for key, vals in feature_mapping.items():
        if isinstance(vals, (list, tuple, set)):
            for v in vals:
                rev[v].append(key)
        elif isinstance(vals, str):
            rev[vals].append(key)
        else:
            # if user accidentally provided nested dict for a key, try to flatten keys
            try:
                for v in vals:
                    rev[v].append(key)
            except Exception:
                pass
    return rev

def gff_to_beds(gff_path, outdir, mapping, single_features=set(), skip_mito=True, overwrite=False):
    os.makedirs(outdir, exist_ok=True)
    revmap = build_reverse_map(mapping)

    # containers: key -> list of bed rows (as tuples)
    beds = defaultdict(list)

    # counters for fallback names to avoid duplicate / empty names
    local_name_counters = defaultdict(int)

    with open(gff_path) as fh:
        for lineno, line in enumerate(fh, start=1):
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                # skip malformed or short lines
                continue
            seqid, src, ftype, start_s, end_s, score, strand, phase, attrs = cols[:9]

            # skip mitochondrial contigs
            if skip_mito and is_mito(seqid):
                continue

            # convert coordinates to ints; GFF is 1-based inclusive
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                # bad coords, skip
                continue
            # BED 0-based half-open
            bed_start = max(0, start - 1)
            bed_end = end

            attrd = parse_attrs(attrs)
            name = attrd.get("Name") or attrd.get("ID") or attrd.get("IDREF") or attrd.get("gene") or ftype
            if not name:
                # fallback name: ftype_counter
                local_name_counters[ftype] += 1
                name = f"{ftype}_{local_name_counters[ftype]}"

            # determine categories to emit to
            cat_keys = revmap.get(ftype, []).copy()
            # direct single feature matches (if present and not in mapping)
            if ftype in single_features and ftype not in revmap:
                cat_keys.append(ftype)

            # if nothing matches, skip (we only emit functional categories)
            if not cat_keys:
                continue

            # ensure score is valid: GFF might have '.'; use 0 if missing
            score_out = score if score != "." else "0"
            # ensure strand is valid
            strand_out = strand if strand in ("+", "-", ".") else "."

            # build bed tuple and append to each matching category
            for key in cat_keys:
                beds[key].append((seqid, bed_start, bed_end, name, score_out, strand_out))

    # write out bed files
    for key, rows in beds.items():
        # sanitize filename
        fname = f"{key}.bed"
        outpath = os.path.join(outdir, fname)
        if os.path.exists(outpath) and not overwrite:
            print(f"Warning: {outpath} exists and overwrite is False. Appending output.", file=sys.stderr)
            mode = "a"
        else:
            mode = "w"
        with open(outpath, mode) as outfh:
            for r in rows:
                outfh.write("\t".join(map(str, r)) + "\n")
        print(f"Wrote {len(rows)} regions to {outpath}")

def main():
    parser = argparse.ArgumentParser(description="Create regionDB BED files from a GFF using a feature mapping.")
    parser.add_argument("gff", 
                       nargs='?', 
                       default="../data/s288c_annotation_genome.gff",
                       help="Input GFF file (default: ../data/s288c_annotation_genome.gff)")
    parser.add_argument("--outdir", "-o", default="regionDB", help="Output directory for BED files (default: regionDB)")
    parser.add_argument("--no-skip-mito", dest="skip_mito", action="store_false", help="Don't skip mitochondrial contigs")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing BED files rather than append")
    args = parser.parse_args()

    # Check if GFF file exists
    if not os.path.exists(args.gff):
        print(f"Error: GFF file not found: {args.gff}")
        print("Please ensure the S. cerevisiae annotation file is available")
        sys.exit(1)
    
    print(f"Using GFF file: {args.gff}")
    print(f"Output directory: {args.outdir}")

    # optionally extend SINGLE_FEATURES from mapping keys that are plain strings (not lists)
    # but we already pre-defined SINGLE_FEATURES above; you can also auto-detect:
    # for k,v in FEATURE_MAPPING.items():
    #     if isinstance(v, str):
    #         SINGLE_FEATURES.add(v)

    gff_to_beds(args.gff, args.outdir, FEATURE_MAPPING, SINGLE_FEATURES, skip_mito=args.skip_mito, overwrite=args.overwrite)

if __name__ == "__main__":
    main()

