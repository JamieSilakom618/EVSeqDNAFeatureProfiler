# EV-seq DNA Feature Profiler

This repository contains a computational workflow for Extracellular Vesicle sequencing (EV-seq) analysis. The workflow enables reproducible analysis of DNA content within extracellular vesicles for *Saccharomyces cerevisiae*.

## Abstract

This computational pipeline analyzes DNA content in extracellular vesicles through a systematic three-step approach: quality control and genome alignment with genomic composition analysis, feature-level enrichment testing, and abundance quantification. The workflow is designed for reproducibility and follows best practices for bioinformatics research.

## workflow Information
- **Author**: Nutticha Silakom
- **Institution**: Chulalongkorn University, Bangkok, Thailand
- **Program**: Bioinformatics and Computational Biology, Graduate School
- **Workflow Version**: 1.0.0
- **Last Updated**: December 2025
- **Data Availability**: All scripts and analysis parameters are provided for full reproducibility

## Workflow Overview

The analysis pipeline consists of three main steps, each corresponding to sections in the manuscript methodology:

1. **Read Processing and Genomic Composition Analysis** (`01_read_processing_and_genomic_composition/`)
2. **Feature-Level Locus Enrichment Analysis** (`02_locus_enrichment/`)
3. **Quantification of EV DNA Abundance** (`03_abundance_quantification/`)

## Requirements

**Minimum System Requirements:**
- Ubuntu 18.04+ / macOS 10.15+ / CentOS 7+
- 16 GB RAM (32 GB recommended)
- 100 GB free disk space

**Software Dependencies:**
- Python 3.8+ (with pip)
- R 4.0+ with Bioconductor
- BWA-MEM 0.7.17+
- SAMtools 1.10+
- BEDtools 2.29.0+
- FastQC 0.11.8+

**Installation:**
```bash
# Install Python dependencies
pip install -r requirements.txt

# Install R dependencies
Rscript -e "install.packages(c('LOLA', 'GenomicRanges'))"
```

## Quick Start

### Input Data Requirements
All input data files should be placed in the `data/` folder:
```
data/
├── gene_fpkm.csv                      # EV-seq FPKM results  
├── srr5658399_count_with_length.csv   # RNA-seq reference data
└── s288c_annotation_genome.gff        # S. cerevisiae annotation
```

### Workflow Steps
1. **Configure Parameters**: Edit `config.ini` with your file paths and parameters
2. **Step 1 - Read Processing**: Process FASTQ files and analyze genomic composition
   ```bash
   cd 01_read_processing_and_genomic_composition/
   ./qc_alignment.sh /path/to/your/fastq/files
   ```
3. **Step 2 - Locus Enrichment**: Perform LOLA enrichment analysis
   ```bash
   cd ../02_locus_enrichment/
   Rscript lola_run.R
   ```
4. **Step 3 - Abundance Quantification**: Calculate FPKM and correlations
   ```bash
   cd ../03_abundance_quantification/
   python calculate_fpkm.py
   python correlation_test.py
   ```

## Directory Structure

```
exo_seq/
├── 01_read_processing_and_genomic_composition/ # Quality control, alignment & composition analysis
├── 02_locus_enrichment/                        # LOLA analysis for genomic feature enrichment
├── 03_abundance_quantification/                # FPKM calculation and correlation analysis
├── data/                                       # Reference genomes and annotations
└── results/                                    # Final output files
```


## References

### Core Tools and Software

1. **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

2. **BWA-MEM**: Li, H. & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754-1760. https://doi.org/10.1093/bioinformatics/btp324

3. **SAMtools**: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics, 25(16), 2078-2080. https://doi.org/10.1093/bioinformatics/btp352

4. **BEDtools**: Quinlan, A. R. & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841-842. https://doi.org/10.1093/bioinformatics/btq033

5. **LOLA**: Sheffield, N. C. & Bock, C. (2016). LOLA: enrichment analysis for genomic region sets and regulatory elements in R and Bioconductor. Bioinformatics, 32(4), 587-589. https://doi.org/10.1093/bioinformatics/btv612


### Statistical Methods

6. **Spearman Correlation**: Spearman, C. (1904). The proof and measurement of association between two things. The American Journal of Psychology, 15(1), 72-101.

7. **Fisher's Exact Test**: Fisher, R. A. (1922). On the interpretation of χ² from contingency tables, and the calculation of P. Journal of the Royal Statistical Society, 85(1), 87-94.

### Genome References

8. **S. cerevisiae Reference Genome**: Cherry, J. M., Hong, E. L., Amundsen, C., Balakrishnan, R., Binkley, G., Chan, E. T., ... & Wong, E. D. (2012). Saccharomyces Genome Database: the genomics resource of budding yeast. Nucleic Acids Research, 40(D1), D700-D705. https://doi.org/10.1093/nar/gkr1029