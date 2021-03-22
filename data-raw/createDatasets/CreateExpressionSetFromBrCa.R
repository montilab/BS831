## generation of the ExpressionSet for Homework1
##
dat <- readRDS(file.path(OMPATH,"data/TCGA.BRCA25.rds"))
dat1 <- variationFilter(dat,ngenes=500,do.plot=TRUE)
pData(dat1) <- pData(dat1)[,c("er_status","pr_status","h2_ihc","h2_fish","h2_status","pam50","subtype")]
write.table(exprs(dat1),file=file.path(OMPATH,"data/tmp/ESet_exprs.xls"),sep="\t",quote=TRUE)
write.table(pData(dat1),file=file.path(OMPATH,"data/tmp/ESet_pData.xls"),sep="\t",quote=TRUE)
write.table(fData(dat1),file=file.path(OMPATH,"data/tmp/ESet_fData.xls"),sep="\t",quote=TRUE)
