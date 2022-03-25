## Creation of different types of datasets for students to "explore"

## settings
library(Biobase)
library(ggplot2)
library(dplyr)
DPATH <- file.path(Sys.getenv("OMPATH"),"data")

##########################################################################
## support functions                                                    ##
##########################################################################
## ANONYMIZE ESET
anonymize_eset <- function(eset,rm_fdata=TRUE, seed=123)
{
    ndigits <- function(x) {floor(log10(x)) + 1}
    sampleNames(eset) <-
        sprintf(paste0("sample_%0",ndigits(ncol(eset)),"d"),sample(ncol(eset)))
    featureNames(eset) <-
        sprintf(paste0("feature_%0",ndigits(nrow(eset)),"d"),sample(nrow(eset)))
    if ( rm_fdata ) {
        fData(eset) <- fData(eset) %>%
            tibble::rownames_to_column() %>%
            dplyr::mutate(feature_name=rowname) %>%
            tibble::column_to_rownames() %>%
            dplyr::select(feature_name)
    }
    return(eset)
}
## IMPUTE MIN VALUES
impute_min_values <- function(dat,shrink=TRUE,seed=1928)
{
    for ( i in 1:nrow(dat) )
    {
        if ( sum( na_idx <- is.na(dat[i,]) )<1 ) next
        upper <- min(dat[i,], na.rm=TRUE)
        lower <- if(shrink) upper/2 else 0.0
        dat[i,na_idx] <- runif(n=sum(na_idx),min=lower,max=upper)
    }
    return(dat)
}
##########################################################################
## A proteomics dataset                                                 ##
##########################################################################
## log(x) is ~normal
eset1 <- readRDS(file.path(Sys.getenv("NECS"),
                           "data_repository/proteomics/somascan_eSet.rds"))
## reduce to 2000 features
eset1 <- eset1[order(matrixStats::rowMads(log2(exprs(eset1))),
                     decreasing = TRUE)[1:1000],]
pData(eset1) <- pData(eset1) %>%
    dplyr::mutate(Cohort=factor(case_when(
        Cohort=="Septuagenarian Control" ~ "Control",
        Cohort=="Spousal Control" ~ "Control",
        Cohort=="Centenarian Offspring" ~ "Offspring",
        Cohort=="Centenarian Singleton" ~ "Centenarian",
        Cohort=="Centenarian Sib-Pair" ~ "Centenarian",
        TRUE ~ ""),levels=c("Control","Offspring","Centenarian")))

## anonymize samples and features
eset1 <- anonymize_eset(eset1,seed=4321)
## additional sample annotation anonymization
set.seed(123)
pData(eset1) <- pData(eset1) %>%
    dplyr::select(Gender,Age,Cohort,smoking,APOE3,outlier) %>%
    dplyr::mutate(Cohort=factor(c("pheno1","pheno3","pheno2")[Cohort])) %>%
    dplyr::mutate(Age=Age-25) %>%
    dplyr::mutate(APOE3=sample(paste0("geno",1:nlevels(APOE3)))[APOE3]) %>%
    dplyr::rename(genotype="APOE3",phenotype="Cohort") %>%
    dplyr::relocate(outlier,.after="genotype")

## check the distribution
par(mfrow=c(1,2))
hist(exprs(eset1),main="raw values")
hist(log2(exprs(eset1)),main="log-transformed values")

## PCA
pca1 <- prcomp(t(exprs(eset1)), scale = TRUE ) ## perform PCA
summary(pca1)$importance[,1:4]                 ## show variance explained by each component
DF1 <- dplyr::inner_join(
    data.frame(pca1$x) %>% tibble::rownames_to_column("sampleID"),
    pData(eset1) %>% dplyr::select(phenotype) %>% tibble::rownames_to_column("sampleID"),
    by="sampleID") %>%
    dplyr::select(sampleID,phenotype,PC1:PC5)
ggplot(DF1,aes(x=PC1,y=PC2,col=phenotype)) +
    geom_point() +
    geom_point(data=DF1 %>% filter(PC1>60),aes(x=PC1,y=PC2),color="red",size=3)

## remove outliers and re-plot
DF1_nooutlier <- DF1 %>% filter(PC1<60)
ggplot(DF1_nooutlier,aes(x=PC1,y=PC2,col=phenotype)) +
    geom_point()

##########################################################################
## A count dataset                                                      ##
##########################################################################
## log(count) is ~normal
## with lots of zero counts

eset2 <- readRDS(file.path(Sys.getenv("OMPATH"),"data/HNSC_htseq_raw_counts_AEvsG1vsG3.RDS"))
par(mfrow=c(1,2))
hist(exprs(eset2))
hist(log2(exprs(eset2)+1))
set.seed(123)
eset2 <- eset2[sample(nrow(eset2),1000),]
par(mfrow=c(1,2))
hist(exprs(eset2))
hist(log2(exprs(eset2)+1))

##########################################################################
## A metabolomics dataset                                               ##
##########################################################################
## log(x) is ~normal
## missing values
eset3 <- readRDS(file.path(Sys.getenv("NECS"),"multiomics/workdir/necs2021metabolon_wmissing.eSet.Rds"))
eset3 <- eset3[order(matrixStats::rowMads(log2(exprs(eset3)),na.rm = TRUE),decreasing = TRUE)[1:1000],]
eset3 <- anonymize_eset(eset3,seed=7890)
## additional sample annotation anonymization
set.seed(123)
pData(eset3) <- pData(eset3) %>%
    dplyr::select(Gender,Age,Cohort) %>%
    dplyr::mutate(Cohort=factor(c("pheno2","pheno3","pheno1")[Cohort])) %>%
    dplyr::mutate(Age=Age-27) %>%
    dplyr::rename(phenotype="Cohort")

## PCA (on imputed data)
eset3_imputed <- eset3
exprs(eset3_imputed) <- impute_min_values(exprs(eset3))

pca3 <- prcomp(t(exprs(eset3_imputed)), scale = TRUE ) ## perform PCA
summary(pca3)$importance[,1:4]                 ## show variance explained by each component
DF3 <- dplyr::inner_join(
    data.frame(pca3$x) %>% tibble::rownames_to_column("sampleID"),
    pData(eset3_imputed) %>% dplyr::select(phenotype) %>% tibble::rownames_to_column("sampleID"),
    by="sampleID") %>%
    dplyr::select(sampleID,phenotype,PC1:PC5)
ggplot(DF3,aes(x=PC1,y=PC2,col=phenotype)) +
    geom_point()
ggplot(DF3,aes(x=PC1,y=PC2,col=phenotype)) +
    geom_point() +
    geom_point(data=DF3 %>% filter(PC1< -35),aes(x=PC1,y=PC2),color="red",size=3)

## remove outliers and re-plot
DF3_nooutlier <- DF3 %>% filter(PC1> -35)
ggplot(DF3_nooutlier,aes(x=PC1,y=PC2,col=phenotype)) +
    geom_point()

##########################################################################
## Save Datasets                                                        ##
##########################################################################

saveRDS(eset1,file=file.path(DPATH,"eset1_lookatthedata.rds"))
saveRDS(eset2,file=file.path(DPATH,"eset2_lookatthedata.rds"))
saveRDS(eset3,file=file.path(DPATH,"eset3_lookatthedata.rds"))


