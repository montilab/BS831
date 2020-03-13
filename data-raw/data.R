library(usethis)

HNSC_RNASeq_ES <- readRDS("data-raw/objects/HNSC_RNASeq_ES.rds")
usethis::use_data(HNSC_RNASeq_ES)

HNSC_RNASeq_toy_ES <- readRDS("data-raw/objects/HNSC_RNASeq_toy_ES.rds")
usethis::use_data(HNSC_RNASeq_toy_ES)

dm10 <- readRDS("data-raw/objects/dm10.rds")
usethis::use_data(dm10)

AEDAT.collapsed.mad4k <- readRDS("data-raw/objects/AEDAT.collapsed.mad4k.rds")
usethis::use_data(AEDAT.collapsed.mad4k)

breast_loi_133p2 <- readRDS("data-raw/objects/breast_loi_133p2.RDS")
usethis::use_data(breast_loi_133p2)

HNSC_htseq_raw_counts_AEvsG1vsG3 <- readRDS("data-raw/objects/HNSC_htseq_raw_counts_AEvsG1vsG3.RDS")
usethis::use_data(HNSC_htseq_raw_counts_AEvsG1vsG3)

zebrafish_cufflinks_counts_fpkm_eSet <- readRDS("data-raw/objects/zebrafish_cufflinks_counts_fpkm_eSet.RDS")
usethis::use_data(zebrafish_cufflinks_counts_fpkm_eSet)

zebrafish_htseq_raw_counts_eSet <- readRDS("data-raw/objects/zebrafish_htseq_raw_counts_eSet.RDS")
usethis::use_data(zebrafish_htseq_raw_counts_eSet)

renamedBreastDB <- readRDS("data-raw/objects/renamedBreastDB.RDS")
usethis::use_data(renamedBreastDB)

AEDAT.collapsed <- readRDS("data-raw/objects/AEDAT.collapsed.RDS")
usethis::use_data(AEDAT.collapsed)

HNSC_htseq_raw_counts_toy_AEvsG2 <- readRDS("data-raw/objects/HNSC_htseq_raw_counts_toy_AEvsG2.rds")
usethis::use_data(HNSC_htseq_raw_counts_toy_AEvsG2)
