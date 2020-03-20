## see scripts at ~/Research/Projects/environcology/DrugMatrix/scripts/
##
require(Biobase)
setwd("~/Research/Projects/environcology/DrugMatrix/")
source("scripts/drugmatrix.support.functions.R")

## LOAD DATASET

data_curated <- readRDS('../datasets/data_curated_mich_ensg.RDS')
DAT <- data_curated[,pData(data_curated)[,"ORGAN.OR.CELL.TYPE"]=="LIVER"]

fname <- "workdir/DAT1liver.RData"
if ( file.access(fname)==0 ) {
    load(file=fname)
}
if ( file.access(fname)!=0 ) {
  DAT1 <- annotate.data(DAT,release="current")
  MAP1 <- annotate.data(DAT,release="current",map=TRUE)
  save(DAT1,MAP1,file=fname)
}
## STEP 2) create a list of list PAIRS (dataset dependent)

pairing <- {
  if (FALSE) # do not collapse multiple doses nor multiple durations
    control.pairing(pData(DAT1),
                    CTLids=c("ORGAN.OR.CELL.TYPE", "DURATION", "ROUTE", "VEHICLE"),
                    TRTids=c("ORGAN.OR.CELL.TYPE", "DURATION", "ROUTE", "VEHICLE", "CHEMICAL", "DOSE"))
  else              # collapse multiple doses and multiple durations
    control.pairing(pData(DAT1),
                    CTLids=c("ORGAN.OR.CELL.TYPE", "ROUTE"),
                    TRTids=c("ORGAN.OR.CELL.TYPE", "ROUTE", "CHEMICAL"),min.ctl=3,min.trt=3)
#    control.pairing(pData(DAT1),
#                    CTLids=c("ORGAN.OR.CELL.TYPE", "ROUTE", "VEHICLE"),
#                    TRTids=c("ORGAN.OR.CELL.TYPE", "ROUTE", "VEHICLE", "CHEMICAL"),min.ctl=3,min.trt=3)
}
## create descriptive names for list pairs
PAIRS1 <- pairing$pairs
sID <- sapply(PAIRS1,function(z) z$treatment[1])
PANN1 <- pData(DAT1)[match(sID,rownames(pData(DAT1))),
                            c("CHEMICAL","DOSE","DURATION","Carcinogen_liv","GenTox")]
tmp <- apply(PANN1[,1:3],1,paste,collapse=' /// ')
if ( length(tmp)!=length(unique(tmp)) ) stop( 'non-unique names' )
names(PAIRS1) <- rownames(PANN1) <- tmp

## extract a subset of 5 non-carcinogen and 5 carcinogen

subset <- PANN1[,"CHEMICAL"] %in% c("ACETAMINOPHEN","ASPIRIN","DOXORUBICIN","IBUPROFEN","ROSIGLITAZONE",                      # non-carcinogens
                                   "2,3,7,8-TETRACHLORODIBENZO-P-DIOXIN","DIAZEPAM","LOVASTATIN","TAMOXIFEN","TESTOSTERONE") # carcinogens

PAIRS2 <- PAIRS1[subset]
PANN2 <- PANN1[subset,]
DAT2 <- DAT1[,unique(unlist(PAIRS2))]
pData(DAT2) <- pData(DAT2)[,c("CHEMICAL","DOSE","DURATION","SEX","Carcinogen_liv","GenTox")]
pData(DAT2)[,"CHEMICAL"] <- factor(pData(DAT2)[,"CHEMICAL"],
                                   levels=levels(pData(DAT2)[,"CHEMICAL"])[levels(pData(DAT2)[,"CHEMICAL"]) %in% unique(pData(DAT2)[,"CHEMICAL"])])
CLS <- pData(DAT2)[,c("Carcinogen_liv")]
pData(DAT2)[,c("Carcinogen_liv")] <- factor(c("NON-CARC","CARCINOGEN")[CLS+1],levels=c("NON-CARC","CARCINOGEN"))
table(pData(DAT2)[,"CHEMICAL"])

## save data

dm10 <- DAT2
dmMap <- PAIRS2
dmAnn <- PANN2
save(dm10,dmMap,dmAnn,file="workdir/DMsubset.RData")

## let's save a small subset properly sorted and filtered

# load(file=file.path(OMPATH,"data/DMsubset.RData"))
dat <- dm10[,pData(dm10)$Carcinogen_liv %in% c('NON-CARC','CARCINOGEN')]
dat <- dat[,order(pData(dat)$Carcinogen_liv)] ## order according to phenotype labels
table(pData(dat)$Carcinogen_liv)              ## show class composition
dat500 <- variationFilter(dat,ngenes=500,score="mad")
saveRDS(dat500,file=file.path(OMPATH,"data/DM10.500genes.RDS"))

##
require(CBMRtools)

DAT3 <- variationFilter(DAT2,ngenes=2000)

hc.row<-hcopt(as.dist(1-cor(t(exprs(DAT3)))),method="ward.D")
hc.col <- hcopt(dist(t(exprs(DAT3))), method="ward.D")

p2 <- heatmap.ggplot2(eSet=DAT3, col.clust = TRUE, row.clust = TRUE,
                      col.clust.hc = hc.col, row.clust.hc = hc.row,
                      col.lab = c("CHEMICAL","Carcinogen_liv","GenTox"), row.lab = "",
                      heatmap.y.text = FALSE, heatmap.x.text = FALSE,
                      heatmap.colorlegend.name = "RNASeq_expression",
                      title.text = "DrugMatrix subset of 10 compounds",
                      col.legend.name = c("CHEMICAL","Carcinogen_liv","GenTox"),
                      row.legend.name = "",
                      row.scaling = "z-score.capped",
                      z.norm = FALSE,
                      cuttree.col = 4, cuttree.row = 3,
                      verbose = FALSE, show = FALSE)
grid.arrange(p2)

