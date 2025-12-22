# Step 2: Feature-Level Locus Enrichment Analysis

This directory contains the LOLA (Locus Overlap Analysis) workflow for assessing genomic feature enrichment.

## Process Overview

1. **Genomic Binning**: EV-DNA reads are merged and partitioned into 150-bp genomic bins
2. **Feature Annotation**: Bins are intersected with annotated features from S. cerevisiae R64-2-1 GFF
3. **Enrichment Testing**: Statistical enrichment assessment using Fisher's exact test
4. **Feature Ranking**: Features ranked by statistical significance and effect size

## Files

- `lola_run.R` - Main LOLA analysis script
- `prepare_mtfeatures_forLOLA_regionDB.py` - Prepare mitochondrial feature database
- `prepare_regionBD_forlola.py` - Prepare region database for LOLA
- `regionDB/` - Region database for different genomic compartments
- `universe/` - Universe sets (binned genomic regions)
- `userset/` - User sets (EV-DNA derived regions)
- `lola_output/` - LOLA results and statistics

## Usage

```bash
# Prepare region databases (uses ../data/s288c_annotation_genome.gff by default)
python prepare_regionBD_forlola.py
python prepare_mtfeatures_forLOLA_regionDB.py

# Or specify custom GFF file
python prepare_regionBD_forlola.py /path/to/custom.gff
python prepare_mtfeatures_forLOLA_regionDB.py /path/to/custom.gff

# Run LOLA analysis
Rscript lola_run.R
```

## Analysis Levels

The analysis is performed at multiple hierarchical levels:
1. Whole-genome (nuclear and mitochondrial)
2. Nuclear genome only
3. Mitochondrial genome only

## Output

- Enrichment statistics for genomic features
- Identification of disproportionately represented elements in EV DNA