require(Biobase)
require(CBMRtools)
OMPATH <- Sys.getenv("OMPATH")

eSet <- readRDS(file.path(OMPATH,"data/AEDAT.collapsed.RDS"))
eSet <- variationFilter(eSet,ngenes=3000,score="mad")
featureNames(eSet) <- fData(eSet)$hgnc_symbol # let's use gene symbols as feature names

eSet1 <- eSet[,pData(eSet)$Characteristics.DiseaseState %in%
               c("non-basal-like breast cancer","sporadic basal-like breast cancer")]
pData(eSet1) <- pData(eSet1)[,c("Characteristics.Individual","Characteristics.DiseaseState")]
colnames(pData(eSet1)) <- c("individual","diseaseState")

eSet1 <- eSet1[,order(pData(eSet1)$diseaseState)]
eSet1$diseaseState <- factor(c("nonbasal","basal")[match(eSet1$diseaseState,unique(eSet1$diseaseState))],
                             levels=c("nonbasal","basal"))
saveRDS(eSet1,file=file.path(OMPATH,"data/renamedBreastDB.RDS"))
    
