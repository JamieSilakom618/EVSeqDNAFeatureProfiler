# Step 1: Read Processing and Genomic Composition Analysis

This directory contains scripts and results for the initial processing of raw sequencing reads and genomic composition analysis.

## Process Overview

### Read Processing
1. **Quality Control**: Raw reads are quality-checked using FastQC
2. **Adapter Trimming**: Adapters are trimmed using standard trimming tools
3. **Genome Alignment**: High-quality reads are aligned to the S288C reference genome using BWA-MEM
4. **Mapping Statistics**: Genome coverage and read distribution assessment using SAMtools

### Genomic Composition Analysis
5. **Read Quantification**: Mapped reads quantified across mitochondrial and nuclear genomes
6. **Compositional Analysis**: Establishes whether EV-DNA pools reflect the genomic composition observed in RNA-seq
7. **Baseline Establishment**: Provides baseline expectation for subsequent analyses

## Files

- `qc_alignment.sh` - Main script for quality control and alignment
- `*.fastqc.html` - FastQC quality control reports
- `*.bam` - Aligned reads (output)
- `*.bai` - BAM index files (output)

## Usage

```bash
# Run quality control and alignment
./qc_alignment.sh input_R1.fastq input_R2.fastq output_prefix

# View alignment statistics
samtools flagstat output_prefix.bam

# Calculate genomic composition
samtools view -c -F 4 output_prefix.bam  # Total mapped reads
samtools view -c -F 4 output_prefix.bam chrMT  # Mitochondrial reads
```

## Requirements

- FastQC
- BWA-MEM
- SAMtools
- S288C reference genome

## Analysis Output

**Genomic Composition Metrics:**
- Proportion of reads mapping to mitochondrial vs nuclear DNA
- Comparison with expected genomic composition from RNA-seq data
- Quality metrics for EV-DNA preparation

## Output

Aligned BAM files and genomic composition statistics are used as input for downstream locus enrichment analysis.