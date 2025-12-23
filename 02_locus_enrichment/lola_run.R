#!/usr/bin/env Rscript
# run_lola_custom.R
# Run LOLA with user-provided userSets, universe, and a custom regionDB.
# Edit the path variables below as needed.

library(LOLA)
library(GenomicRanges)
library(utils)
# -----------------------------
# User-editable paths
# -----------------------------

wg_user  = readBed('./userset/userset_wg_binned.bed')
wg_uni = readBed('./universe/binned150_allgenome.bed')
regionDB_wg <- loadRegionDB("regionDB", collections = "wg")
res_wg <- runLOLA(wg_user, wg_uni, regionDB_wg)
checkUniverseAppropriateness(wg_user, wg_uni)


regionDB_nuc <- loadRegionDB("regionDB", collections = "chr_feature")
res_nuc <- runLOLA(wg_user, wg_uni, regionDB_nuc)


regionDB_mt <- loadRegionDB("regionDB", collections = "mt_feature")
res_mt <- runLOLA(wg_user, wg_uni, regionDB_mt)


#head(res_wg)

write.csv(res_wg, file = "LOLA_genomic_results_validated.csv", row.names = FALSE)
write.csv(res_nuc, file = "LOLA_chr_results_validated.csv", row.names = FALSE)
write.csv(res_mt, file = "LOLA_mt_results_validated.csv", row.names = FALSE)
