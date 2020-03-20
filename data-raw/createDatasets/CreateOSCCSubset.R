require(Biobase)

setwd("~/Research/Presentations/2015_10_09_BS830/")

OSCC <- readRDS("~/Research/Projects/oralcancer/tcga/firehose_2014_12_06/TCGA_OSCC_mRNA.annotated.rds")


colSums(is.na(pData(OSCC))[,c("patient.follow_ups.follow_up-2.days_to_last_followup",
                              "patient.follow_ups.follow_up-3.days_to_last_followup",
                              "patient.follow_ups.follow_up-4.days_to_last_followup")])
TMP <- pData(OSCC)[,c("patient.follow_ups.follow_up-2.days_to_last_followup",
                      "patient.follow_ups.follow_up-3.days_to_last_followup",
                      "patient.follow_ups.follow_up-4.days_to_last_followup")]
TMP[,] <- apply(TMP,2,as.numeric)
colnames(TMP) <- 2:4

FU <- apply((TMP),1,max,na.rm=TRUE)
cbind(TMP,FU)[1:100,]
