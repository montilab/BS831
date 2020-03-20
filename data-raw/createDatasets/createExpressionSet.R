## Generating a simple ExpressionSet for the students to experiment with

require(Biobase)
OMPATH <- Sys.getenv("OMPATH"); print(OMPATH)

set.seed(123)
nS <- 10
nG <- 100

EXP <- matrix(rlnorm(nS*nG),nrow=nG,ncol=nS,
              dimnames=list(sprintf("gene%03d",1:nG),
                  sprintf("sample%02d",1:nS)))

FDAT <- data.frame(rownames(EXP),
                   gDescription1=sprintf("gDescription1.%03d",1:nG),
                   gDescription2=sprintf("gDescription2.%03d",1:nG),
                   row.names=1)

PDAT <- data.frame(colnames(EXP),
                   sDescription1=sprintf("sDescription1.%02d",1:nS),
                   sDescription2=sprintf("sDescription2.%02d",1:nS),
                   sDescription3=sprintf("sDescription3.%02d",1:nS),
                   row.names=1)

ESET <- ExpressionSet(assayData=EXP,
                      phenoData=AnnotatedDataFrame(PDAT),
                      featureData=AnnotatedDataFrame(FDAT))

## scramble EXP row and column order
EXP1 <- EXP[sample(1:nrow(EXP)),]

write.table(EXP1,file=file.path(OMPATH,"data/ESet_exprs.xls"))
write.table(FDAT,file=file.path(OMPATH,"data/ESet_fData.xls"))
write.table(PDAT,file=file.path(OMPATH,"data/ESet_pData.xls"))

