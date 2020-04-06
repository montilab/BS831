require(Biobase)
source( file.path(Sys.getenv("CBMRtools"),"CBMRtools/R/variationFilter.R") )
source( file.path(Sys.getenv("CBMRtools"),"CBMRtools/R/misc.R") )
OMPATH <- Sys.getenv("OMPATH")

## researchData is equivalent to CBMrepositoryData on SCC
brca <- readRDS("~/researchData/tcga/esets_processed/BRCA/BRCA_2018_11_01_subtyped_DESeq2_eSet.rds")
brca$pam50 <- factor(brca$pam50,levels=c("LumA","LumB","Her2","Basal","Normal"))
exprs(brca) <- log2(exprs(brca)+1)

## extract 25 samples from each subtype (excluding ("Normal")
N <- 25
set.seed(123)
sIDX <- c(sample(sampleNames(brca[,brca$pam50=="LumA"]),size=N,replace=FALSE),
          sample(sampleNames(brca[,brca$pam50=="LumB"]),size=N,replace=FALSE),
          sample(sampleNames(brca[,brca$pam50=="Her2"]),size=N,replace=FALSE),
          sample(sampleNames(brca[,brca$pam50=="Basal"]),size=N,replace=FALSE))

brca25  <- variationFilter(brca,ngenes=5000,score="mad",do.plot=TRUE)[,sIDX]
brca25$pam50 <- droplevels(brca25$pam50)
saveRDS( brca25, file.path(OMPATH,"data/TCGA.BRCA25.rds") )
save(brca25,"data/TCGA.BRCA25.rda")
