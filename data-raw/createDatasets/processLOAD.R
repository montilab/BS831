## Read in LOAD data
LOAD <- readRDS( file.path(OMPATH,"data/LOAD1.RDS"))

## restrict to probes with mapping gene symbols
LOAD1 <- LOAD[!is.na(fData(LOAD)$gene_symbol) & fData(LOAD)$gene_symbol!="",]; nrow(LOAD1)

## let us make the phenotype annotation more intuitive
if(FALSE) {
pData(LOAD1) <- pData(LOAD1) %>%
    dplyr::rename(loadStatus="characteristics_ch2") %>%
    dplyr::rename(brainTissue="characteristics_ch2.8") %>%
    dplyr::mutate(loadStatus=as.character(loadStatus)) %>%
    dplyr::select(loadStatus,brainTissue) %>%
    dplyr::mutate(newVariable = dplyr::case_when(
        loadStatus=="disease: A" ~ "load",
        loadStatus=="disease: N" ~ "control",
        TRUE ~ stop("unrecognized loadStatus:",loadStatus))
    )
}
colnames(pData(LOAD1)) <- gsub("characteristics_ch2","loadStatus",colnames(pData(LOAD1)))
LOAD1$loadStatus <- factor(c("control","load")[as.numeric(LOAD1$loadStatus=="disease: A")+1],
                           levels=c("load","control"))
colnames(pData(LOAD1)) <- gsub("description","brainTissue",colnames(pData(LOAD1)))
print(table(LOAD1$loadStatus,LOAD1$brainTissue))

## some gene symbols have multiple entries
c(total=nrow(LOAD1),unique=length(unique(fData(LOAD1)$gene_symbol)))

## we 'uniquefy' by taking the replicate w/ largest variation
## (one can also use dplyr to do the same perhaps more elegantly)
SD <- apply(exprs(LOAD1),1,sd,na.rm=TRUE)
LOAD1 <- LOAD1[order(SD,decreasing=TRUE),]
LOAD1 <- LOAD1[match(unique(fData(LOAD1)$gene_symbol),fData(LOAD1)$gene_symbol),]
featureNames(LOAD1) <- fData(LOAD1)$gene_symbol
c(total=nrow(LOAD1),unique=length(unique(fData(LOAD1)$gene_symbol)))

## remove probes w/ NA's (alternative would be to impute them)
LOAD2 <- LOAD1[!apply(is.na(exprs(LOAD1)),1,any),]
c(missing.LOAD1=sum(is.na(exprs(LOAD1))),missing.LOAD2=sum(is.na(exprs(LOAD2))))

## save cleaned data
saveRDS(LOAD2,file=file.path(OMPATH,"data/LOAD2.RDS"))
