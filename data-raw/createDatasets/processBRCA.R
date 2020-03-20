## Clean up of raw counts and CPM normalized data from HTset of the HNSC data.
## This is done in order to simplify the pData, as well as to add the fData to the ExpressionSet's

require(Biobase)
require(biomaRt)
require(CBMRtools)

source(paste(Sys.getenv("CBMDEV"),"R/probesetAnnotation.R",sep="/"))
OMPATH <- Sys.getenv("OMPATH")

setwd(Sys.getenv("OMPATH"))

## load RPKM normalized data
brca <- readRDS("~/researchData/tcga/ESets/RNAseq/BRCA_2016_01_28_ES.wPAM50.rds")
pData(brca) <- pData(brca)[,c("subtype",
                              "er_status",
                              "pr_status",
                              "her2_status",
                              "her2_ihc",
                              "her2_fish")]
#                              "patient.age_at_initial_pathologic_diagnosis",
#                              "patient.anatomic_neoplasm_subdivision",
#                              "patient.clinical_cqcf.tumor_type",
#                              "patient.ethnicity",
#                              "patient.follow_ups.follow_up.vital_status",
#                              "patient.follow_ups.follow_up.year_of_form_completion",
#                              "patient.frequency_of_alcohol_consumption",
#                              "patient.gender",
#                              "patient.histological_type",
#                              "patient.hpv_status_by_ish_testing",
#                              "patient.hpv_status_by_p16_testing",
#                              "patient.lymph_node_examined_count",
#                              "patient.neoplasm_histologic_grade",
#                              "patient.number_of_lymphnodes_positive_by_he",
#                              "patient.number_of_lymphnodes_positive_by_ihc",
#                              "patient.number_pack_years_smoked",
#                              "patient.p53_gene_analysis",
#                              "patient.race",
#                              "patient.stage_event.clinical_stage",
#                              "patient.stopped_smoking_year",
#                              "patient.tumor_tissue_site",
#                              "patient.vital_status",
#                              "patient.year_of_form_completion",
#                              "patient.year_of_initial_pathologic_diagnosis",
#                              "patient.year_of_tobacco_smoking_onset")]
table(pData(brca)$subtype)
brca$subtype <- factor(brca$subtype,levels=sort(unique(brca$subtype)))

## anonymize them
set.seed(123)
NEWidx <- sample(1:ncol(brca),ncol(brca))

brca$orig <- sampleNames(brca)
brca1 <- brca[,NEWidx]
sampleNames(brca1) <- sprintf("sample%04s",1:ncol(brca))
MAP <- data.frame(sampleNames(brca1),brca1$orig)
pData(brca1) <- pData(brca1)[,-match("orig",colnames(pData(brca1)))]

## reorder and subset (to reduce sample size)
brca2 <- brca1
brca2 <- brca1[,order(brca2$subtype)]
set.seed(123)
subIDX <- unlist(sapply(levels(brca2$subtype),function(Z)
    sample(sampleNames(brca2)[brca2$subtype==Z],size=min(50,sum(brca2$subtype==Z)))))
brca2 <- brca2[,subIDX]

saveRDS(brca2,file="data/BRCA_2016_01_28_ES.wPAM50.anonymized.RDS")



