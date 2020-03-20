## Clean up of raw counts and CPM normalized data from HTset of the HNSC data.
## This is done in order to simplify the pData, as well as to add the fData to the ExpressionSet's

require(Biobase)
require(biomaRt)
require(CBMRtools)

source(paste(Sys.getenv("CBMDEV"),"R/probesetAnnotation.R",sep="/"))
OMPATH <- Sys.getenv("OMPATH")

setwd(Sys.getenv("OMPATH"))

## process HTseq raw
raw <- readRDS(file.path(OMPATH,"data/HNSC_htseq_raw_counts.RDS"))
ref <- readRDS("~/Research/CBMrepositoryData/tcga/ESets/HNSC_2015_11_01_ES.rds")

shortenName <- function( X ) { gsub("A$","",paste(unlist(strsplit( X, "\\_" ))[1:4],collapse="-")) }
sampleNames(raw) <- sapply( sampleNames(raw), shortenName )
length(cmn <- intersect(sampleNames(raw),sampleNames(ref)))

table(pData(raw)$patient.stage_event.clinical_stage,useNA="ifany")
table(pData(raw)$patient.neoplasm_histologic_grade,useNA="ifany")

## remove NAs and "gx"
idx <- !is.na(pData(raw)$patient.stage_event.clinical_stage) &
       !is.na(pData(raw)$patient.neoplasm_histologic_grade) &
       pData(raw)$patient.neoplasm_histologic_grade!="gx"
raw1 <- raw[,idx]

## identify the Adjagcent Epithelial (AEs)
ANidx <- grep("-11A-",pData(raw1)$sample_name)
pData(raw1)$patient.stage_event.clinical_stage[ANidx] <- "AE"
pData(raw1)$patient.neoplasm_histologic_grade[ANidx] <- "AE"

## consolidate stage iv
S4idx <- grep("stage iv",pData(raw1)$patient.stage_event.clinical_stage)
pData(raw1)$patient.stage_event.clinical_stage[S4idx] <- "stage iv"

## simplify pData
pdat <- pData(raw1)[,c("patient.age_at_initial_pathologic_diagnosis",
                       "patient.anatomic_neoplasm_subdivision",
                       "patient.clinical_cqcf.tumor_type",
                       "patient.ethnicity"                                                    ,
                       "patient.follow_ups.follow_up.vital_status",
                       "patient.follow_ups.follow_up.year_of_form_completion",
                       "patient.frequency_of_alcohol_consumption",
                       "patient.gender",
                       "patient.histological_type",
                       "patient.hpv_status_by_ish_testing",
                       "patient.hpv_status_by_p16_testing",
                       "patient.lymph_node_examined_count",
                       "patient.neoplasm_histologic_grade",
                       "patient.number_of_lymphnodes_positive_by_he",
                       "patient.number_of_lymphnodes_positive_by_ihc",
                       "patient.number_pack_years_smoked",
                       "patient.p53_gene_analysis",
                       "patient.race",
                       "patient.stage_event.clinical_stage",
                       "patient.stopped_smoking_year",
                       "patient.tumor_tissue_site",
                       "patient.vital_status",
                       "patient.year_of_form_completion",
                       "patient.year_of_initial_pathologic_diagnosis",
                       "patient.year_of_tobacco_smoking_onset")]
pData(raw1) <- pdat

## turn into factors

raw1$patient.neoplasm_histologic_grade) <- as.factor(raw1$patient.neoplasm_histologic_grade)
raw1$patient.stage_event.clinical_stage <- as.factor(raw1$patient.stage_event.clinical_stage)

## scramble order (to increase 'anonymity')
set.seed(123)
NEWidx <- sample(1:ncol(raw1),ncol(raw1))

raw2 <- raw1[,NEWidx]
MAP <- data.frame(sampleNames(raw3),sprintf("sample%03s",1:ncol(raw3)))
sampleNames(raw2) <- MAP[,2]

raw3 <- probesetAnnotation(eset=raw2)

saveRDS(raw3,file="data/HNSC_htseq_raw_counts.anonymized.RDS")
saveRDS(MAP,file="data/HNSC_htseq_anonymizedMAP.RDS")

## pick 1st instance of each ensemble occurrence (duplicates are due
## to multiple entrez IDs matching the same ensemble ID)
##
GS <- GS[match(unique(GS[,"ensembl_gene_id"]),GS[,"ensembl_gene_id"]),]
rownames(GS) <- GS[,'ensembl_gene_id']

## repeat for normalized

nrm <- readRDS("data/HNSC_htseq_normalized_counts.RDS")
nrm1 <- nrm[,sampleNames(raw1)]
all(rownames(pData(nrm1))==rownames(pData(raw1)))
pData(nrm1) <- pData(raw1)

nrm2 <- nrm1[,NEWidx]
all(sampleNames(nrm2)==MAP[,1])
sampleNames(nrm2) <- MAP[,2]

nrm3 <- probesetAnnotation(eset=nrm2)

saveRDS(nrm3,file="data/HNSC_htseq_normalized_counts.anonymized.RDS")

## compare raw to normalized

rawMN <- rowMeans(exprs(raw2))
nrmMN <- rowMeans(exprs(nrm2))

plot(rawMN,nrmMN,pch=".")

################################
## create a simpler (toy) subset
################################

DAT <- readRDS(file=file.path(OMPATH,"data/HNSC_htseq_raw_counts.anonymized.RDS"))
CHM <- c(1:22,"X","Y")

## genes to keep
gIDX <- fData(DAT)[,"chromosome_name"] %in% CHM

## samples to keep (30 samples properly stratified)
set.seed(123)
sIDX <- c(sample(which(DAT$patient.neoplasm_histologic_grade %in% "AE"),15),
          sample(which(DAT$patient.neoplasm_histologic_grade %in% "g2"),15))

TOY <- DAT[gIDX,sIDX]
table(TOY$patient.neoplasm_histologic_grade)

saveRDS(TOY,file=file.path(OMPATH,"data/HNSC_htseq_raw_counts_toy_AEvsG2.rds"))

###########################################
## create a simpler grade1 vs grade3 subset
###########################################

RAW <- readRDS(file=file.path(OMPATH,"data/HNSC_htseq_raw_counts.anonymized.RDS"))
CPM <- readRDS( file.path(OMPATH,"data/HNSC_htseq_normalized_counts.anonymized.RDS") )
CHM <- c(1:22,"X","Y")

## simplify pData labels
colnames(pData(RAW)) <- gsub("patient.neoplasm_histologic_grade","grade",colnames(pData(RAW)))
colnames(pData(RAW)) <- gsub("patient.stage_event.clinical_stage","stage",colnames(pData(RAW)))

## genes to keep
gIDX <-
    ( fData(RAW)[,"chromosome_name"] %in% CHM ) &
        ( fData(RAW)[,"hgnc_symbol"]!="" ) &
            ( !is.na(fData(RAW)[,"hgnc_symbol"]) )


## samples to keep (30 samples properly stratified)
set.seed(123)
sIDX <- c(sample(which(RAW$grade %in% "AE"),40),
          sample(which(RAW$grade %in% "g1"),40),
          sample(which(RAW$grade %in% "g3"),40))

RAW120 <- RAW[gIDX,sIDX]
RAW120$grade <- droplevels(RAW120$grade)
table(RAW120$grade)

CPM120 <- CPM[gIDX,sIDX]
colnames(pData(CPM120)) <- colnames(pData(RAW120))
CPM120$grade <- droplevels(CPM120$grade)
    
saveRDS(RAW120,file=file.path(OMPATH,"data/HNSC_htseq_raw_counts_AEvsG1vsG3.RDS"))
saveRDS(CPM120,file=file.path(OMPATH,"data/HNSC_htseq_normalized_AEvsG1vsG3.RDS"))
