#!/usr/bin/env Rscript
# LOLA Enrichment Analysis for EV-seq Data
#
# This script performs Locus Overlap Analysis (LOLA) to identify genomic features
# enriched in extracellular vesicle DNA.
#
# Authors: Nutticha Silakom
# Institution: Chulalongkorn University, Bangkok, Thailand
# Program: Bioinformatics and Computational Biology, Graduate School
# Version: 1.0.0
# Date: December 2025
#
# Analysis Overview:
# - Performs enrichment testing across multiple genomic compartments
# - Uses Fisher's exact test for statistical significance
# - Analyzes whole-genome, nuclear, and mitochondrial regions separately
#
# Requirements:
# - R >= 4.0.0
# - LOLA >= 1.22.0
# - GenomicRanges >= 1.44.0
#
# Input:
# - User sets: EV-DNA derived genomic regions (userset/)
# - Universe sets: Background genomic regions (universe/)
# - Region databases: Annotated genomic features (regionDB/)
#
# Output:
# - Enrichment statistics for each genomic feature class
# - Statistical significance values and effect sizes

# Load required libraries
suppressPackageStartupMessages({
    library(LOLA)
    library(GenomicRanges)
    library(utils)
})

cat("=====================================\n")
cat("LOLA Enrichment Analysis for EV-seq\n")
cat("Journal of Extracellular Vesicles - Supporting Code\n")
cat("Version 1.0.0\n")
cat(paste("Execution time:", Sys.time(), "\n"))
cat("=====================================\n")

# Set reproducible seed for any random operations
set.seed(12345)

# -----------------------------
# Analysis 1: Whole-genome analysis
# -----------------------------
cat("Running whole-genome analysis...\n")

wg_user  = readBed('./userset/userset_wg_binned.bed')
wg_uni = readBed('./universe/binned150_allgenome.bed')
regionDB_wg <- loadRegionDB("regionDB", collections = "wg")
res_wg <- runLOLA(wg_user, wg_uni, regionDB_wg)
 

nuc_user  = readBed('./userset/userset_nuc.bed')
nuc_uni = readBed('./universe/binned150_nuc.bed')
regionDB_nuc <- loadRegionDB("regionDB", collections = "chr_feature")
res_nuc <- runLOLA(nuc_user, nuc_uni, regionDB_nuc)

mt_user  = readBed('./userset/userset_mt.bed')
mt_uni = readBed('./universe/binned150_mt.bed')
regionDB_mt <- loadRegionDB("regionDB", collections = "mt_feature")
res_mt <- runLOLA(mt_user, mt_uni, regionDB_mt)


#head(res_wg)

write.csv(res_wg, file = "LOLA_genomic_results.csv", row.names = FALSE)
write.csv(res_nuc, file = "LOLA_chr_results.csv", row.names = FALSE)
write.csv(res_mt, file = "LOLA_mt_results.csv", row.names = FALSE)
