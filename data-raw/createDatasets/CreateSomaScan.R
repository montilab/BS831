## CreateSomaScan.R

soma1 <- soma <- readRDS(
    file.path(Sys.getenv("MLAB"),
              "projects/longevity/novartis/data/somascan/180308_CentCollabBU_HybNormPlateScaleMedNormCal.eset.RDS"))

## simplify phenotype
cohort <- as.character(pData(soma)$Cohort)
cohort[grep("Control",cohort)] <- "Control"
cohort <- gsub("Centenarian Offspring","Offspring",cohort)
cohort <- gsub("Centenarian Singleton","Centenarian",cohort)
cohort <- gsub("Centenarian Sib-Pair","Centenarian",cohort)
soma1$Cohort <- factor(cohort,levels=c("Control","Offspring","Centenarian"))
print(table(soma1$Cohort))

## de-identify samples and proteins
soma1 <- soma1[match(unique(fData(soma1)[,"swissprot"]),fData(soma1)[,"swissprot"]),]
featureNames(soma1) <- sprintf("protein%04d",seq(1,nrow(soma1)))
sampleNames(soma1) <- sprintf("sample%03d",seq(1,ncol(soma1)))
fData(soma1) <- fData(soma1)[,1,drop=FALSE]
fData(soma1)[,1] <- featureNames(soma1)
pData(soma1) <- pData(soma1)[,c("Gender","Cohort")]
print(table(pData(soma1)))

soma2 <- soma1[,order(soma1$Cohort)]

## save final object
saveRDS( soma2, file.path(Sys.getenv("OMPATH"),"data/somascan.RDS") )

## show within-sample distribution (subsample)
boxplot(exprs(soma2)[sample(seq(1,nrow(soma2)),size=100),],log="y",las=2,pch="-",
        col=c("gray","pink","red")[soma2$Cohort])
