# EV-seq DNA Feature Profiler

This repository contains a computational workflow for Extracellular Vesicle sequencing (EV-seq) analysis. The workflow enables reproducible analysis of DNA content within extracellular vesicles for *Saccharomyces cerevisiae*.

## Overview

This computational pipeline analyzes DNA content in extracellular vesicles through a systematic three-step approach: quality control and genome alignment with genomic composition analysis, feature-level enrichment testing, and abundance quantification. The workflow is designed for reproducibility and follows best practices for bioinformatics research.

## Workflow Information
- **Author**: Nutticha Silakom
- **Institution**: Chulalongkorn University, Bangkok, Thailand
- **Program**: Bioinformatics and Computational Biology, Graduate School
- **Workflow Version**: 1.0.0
- **Last Updated**: June 2026
- **Data Availability**: All scripts and analysis parameters are provided for full reproducibility
- **License**: MIT License

## Data and Resource Availability

Raw FASTQ and aligned BAM files are not included in this repository due to their large file size. The EV-associated DNA sequencing data generated in this study have been deposited in the DNA Data Bank of Japan (DDBJ) Sequence Read Archive under **BioProject accession PRJDB40093**. The raw sequencing data are available under **Run accession DRR1066343**.

All scripts, source code, analysis workflows, and documentation required to reproduce the analyses presented in this study are provided in this repository.


Once raw data is available:
1. Align reads using the script in `01_read_processing_and_genomic_composition/`:
   ```bash
   cd 01_read_processing_and_genomic_composition/
   bash qc_alignment.sh <sample_id> <fastq_r1> <fastq_r2> <reference_fasta> [threads]
   ```
2. Place the output BAM in the project root:
   ```bash
   cp results/alignment_out/<sample_id>.sorted.bam ../aligned.mapped.sorted.bam
   ```

All downstream output files (FPKM tables, CPM tables, top-1% gene lists) are already included under `03_abundance_quantification/` for reproducibility.

## Workflow Overview

The analysis pipeline has three stages:

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
Place shared input data files in the `data/` folder:
```
data/
├── gene_fpkm.csv                       # EV-seq FPKM results
├── srr5658399_count_with_length.xlsx   # RNA-seq reference count data
└── s288c_annotation_genome.gff         # Reference genome annotation
```

### Workflow Steps
1. **Configure Parameters**: Edit `config.ini` only if you need custom paths
2. **Step 1 - Read Processing**: Run QC, alignment, and BAM generation
   ```bash
   cd 01_read_processing_and_genomic_composition/
   bash qc_alignment.sh <sample_id> <fastq_r1> <fastq_r2> <reference_fasta> [threads]
   ```
   This runs FastQC, aligns with BWA-MEM, and writes `results/alignment_out/<sample_id>.sorted.bam`.
3. **Step 2 - Locus Enrichment**: Prepare region databases and run LOLA
   ```bash
   cd ../02_locus_enrichment/
   python3 prepare_regionBD_forlola.py
   python3 prepare_mtfeatures_forLOLA_regionDB.py
   Rscript lola_run.R
   ```
4. **Step 3 - Abundance Quantification**: Merge EV/SRR data, compute top 1% lists, and annotate genes
   > **Note**: `bash script/run_pipeline.sh` skips steps [1/6] and [2/6] when `aligned.mapped.sorted.bam` is missing, then continues with the Python steps using the outputs already included in the repository.
   ```bash
   cd ../03_abundance_quantification/
   # Run full pipeline (recommended)
   bash script/run_pipeline.sh

   # Or run tasks individually:
   bash task_01_fpkm/script/run_bedtools_coverage.sh
   python3 task_01_fpkm/script/calculate_fpkm.py
   python3 task_02_high_expression/script/merge_file.py          # generates merged CPM + top 1% files
   python3 task_02_high_expression/script/correlation_test.py
   python3 task_02_high_expression/script/intersect_genes.py     # uses default paths; or pass explicit paths
   python3 task_03_gene_annotation/script/sgd_lookup.py         # default input is gene_mt_fpkm.csv
   ```
   
## Key Outputs
The workflow generates the following outputs:

| Stage | Output |
|---|---|
| Read Processing | Sorted BAM files and alignment statistics |
| Genomic Composition | Feature composition summaries |
| Locus Enrichment | LOLA enrichment results |
| Abundance Quantification | Feature-level FPKM and CPM tables |
| High-Abundance Analysis | Top 1% EV and RNA reference gene lists |
| Annotation | SGD gene name annotations |

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