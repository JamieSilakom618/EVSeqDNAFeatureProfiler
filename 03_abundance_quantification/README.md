# Step 3: Quantification of EV DNA Abundance

This directory contains scripts for quantifying genomic feature abundance and correlation analysis.

## Process Overview

1. **Read Counting**: BEDtools coverage-based read quantification normalized to feature length
2. **FPKM Calculation**: Fragments Per Kilobase Million calculation for abundance estimation
3. **High-Expression Analysis**: Identification and analysis of top 1% abundant genes
4. **Correlation Analysis**: Spearman correlation between EV-seq and RNA-seq abundance
5. **Gene Intersection**: Comparative analysis of highly abundant genes across datasets
6. **Statistical Testing**: Assessment at multiple hierarchical genomic levels

## Files

### FPKM Analysis
- `calculate_fpkm.py` - Main FPKM calculation script
- `run_bedtools_coverage.sh` - BEDtools coverage analysis
- `s288c_annotation_genome.gff` - Genome annotation file
- `coverage_count/` - Raw count files for different genomic features
- `fpkm_out/` - FPKM results for all features
- `region_for_fpkm/` - BED files defining genomic regions
- `geneID_search/` - Gene ID lookup and annotation

### High-Expression Gene Analysis
- `compare_gene_high_express/` - Analysis of highly abundant genes
  - `correlation_test.py` - Spearman correlation analysis
  - `merge_file.py` - Merge EV-seq and RNA-seq data
  - `merged_CPM_table.csv` - Merged CPM data for analysis
  - `intersect_genes.py` - Gene intersection analysis script
  - `gene_99thpercentile_high_abundance.csv` - Top 1% highly abundant genes
  - `gene_srr_99thpercetle_high_abundance.csv` - Top 1% genes from SRR data
  - `intersect_genes.csv` - Gene intersection results
  - `intersect_genes_with_names.csv` - Annotated gene intersection results

## Usage

```bash
# Calculate coverage and FPKM
./run_bedtools_coverage.sh
python calculate_fpkm.py

# Perform high-expression gene analysis
cd compare_gene_high_express/
python merge_file.py
python correlation_test.py
python intersect_genes.py
```

## Analysis Levels

Quantification is performed at four hierarchical levels:
1. Whole-genome (nuclear and mitochondrial)
2. Nuclear genome
3. Mitochondrial genome  
4. Individual genomic features

## Features Analyzed

- Genes (protein-coding)
- Non-coding RNAs (ncRNA, rRNA, tRNA, snoRNA, snRNA)
- Regulatory elements (centromeres, telomeres, replication origins)
- Transposable elements
- Pseudogenes
- Mating loci

## Output

- FPKM values for all genomic features
- Correlation coefficients between EV-seq and RNA-seq
- Top 1% highly abundant gene lists
- Gene intersection analysis results
- Feature abundance rankings
- Statistical significance testing results