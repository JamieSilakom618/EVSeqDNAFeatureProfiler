# Changelog

All notable changes to the EV-seq analysis workflow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-12-21

### Added
- Initial release of EV-seq analysis workflow
- Four-step analysis pipeline matching published methodology:
  1. Read Processing and Genome Alignment
  2. Genomic Composition Analysis  
  3. Feature-Level Locus Enrichment Analysis
  4. Quantification of EV DNA Abundance
- Comprehensive documentation for each analysis step
- Master workflow script (`master_workflow.sh`)
- Configuration file (`config.ini`) for centralized parameters
- Requirements file (`requirements.txt`) for Python dependencies
- Citation documentation (`CITATIONS.md`)
- Quality control and alignment scripts
- FPKM calculation and correlation analysis tools
- LOLA enrichment analysis pipeline
- Genomic feature annotation and region databases

### Infrastructure
- Organized directory structure by analysis steps
- README files for each workflow component
- .gitignore file for version control
- Executable workflow scripts
- Python configuration integration

### Dependencies
- Python 3.8+ with scientific libraries (pandas, scipy, numpy)
- R with LOLA package for enrichment analysis
- BWA-MEM for sequence alignment
- SAMtools for BAM file processing
- BEDtools for genomic interval operations
- HTSeq for read quantification
- FastQC for quality control

## [Unreleased]

### Planned
- Add genomic composition analysis scripts for step 2
- Integration with container systems (Docker/Singularity)
- Automated testing framework
- Performance optimizations
- Additional statistical analysis options
- Interactive visualization components