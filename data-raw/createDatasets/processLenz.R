require(biomaRt)
PATH <- "~/dbox/pavia2016/"

load(file="~/bucomp/SCCnb/projects/lymphoma/LymphoClass/datasets/lymphoma_lenz.frma.ensg.RData")
if ( !exists('eSet') ) stop( "!exists('eSet')" )

## annotate probes

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(filters="ensembl_gene_id",
                 attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol","chromosome_name",
                     "band","strand","start_position","end_position","description"),
                 values=featureNames(eSet),
                 mart= mart)
head(mapping)

rMatch <- match(featureNames(eSet), mapping$ensembl_gene_id)
fdat <- matrix(NA,nrow=nrow(eSet),ncol=ncol(mapping),dimnames=list(featureNames(eSet),colnames(mapping)))
fdat[!is.na(rMatch),] <- as.matrix(mapping)[rMatch[!is.na(rMatch)],]
fdat <- data.frame(fdat,stringsAsFactors=FALSE)
fData(eSet) <- fdat

saveRDS(eSet,file=paste(PATH,"data/lymphoma_lenz.frma.ensg.RDS",sep=""))
