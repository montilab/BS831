require(Biobase)    # for ExpressionSet data objects
require(biomaRt)    # for gene annotation
require(CBMRtools)  # 
OMPATH <- Sys.getenv("OMPATH")

## RAW COUNTS (for DEseq2 and edgeR)

## loading htseq raw counts matrix for DEseq2/edgeR analysis
htseq_counts <- read.table(file=file.path(OMPATH,"data/zebrafish_htseq_raw_counts.txt"),
                           htseq_counts_file, sep = "\t", header=TRUE, row.names=1)

## get sample annotation (ignore details and take it as given)
sample_info <- read.table(file=file.path(OMPATH,"data/zebrafish_sample_annotation.txt"),
                          sep="\t", header=TRUE)
sample_info$sample_name_corrected <- gsub("-", "\\.", sample_info$sample_name)
sample_info_order <- match(colnames(htseq_counts), sample_info$sample_name_corrected)
sample_info <- sample_info[sample_info_order,]
rownames(sample_info) <- sample_info$sample_name_corrected

## zebrafish ensembl (Danio rerio is the formal name of https://en.wikipedia.org/wiki/Zebrafish)
mart <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="drerio_gene_ensembl")
mart_map <- getBM(attributes=c("ensembl_gene_id","external_gene_name","hgnc_symbol","wikigene_description"),
                  filter=("ensembl_gene_id"),
                  values = rownames(htseq_counts),
                  mart=mart)

## create ExpressionSet from raw counts (for use w/ DESeq and edgeR)
match_ids <- match(rownames(htseq_counts), mart_map[, "ensembl_gene_id"]) # if multiple matches, get first
htseq_counts1 <- htseq_counts[!is.na(match_ids),]
match_ids1 <- match_ids[!is.na(match_ids)]
if ( any(rownames(htseq_counts1)!=mart_map[match_ids1,"ensembl_gene_id"]) )
    stop( "rownames(htseq_counts1)!=mart_map[match_ids1,\"ensembl_gene_id\"]" )

fdat <- cbind(ensembl_id1=rownames(htseq_counts1),
              mart_map[match_ids1,c("external_gene_name", "hgnc_symbol","wikigene_description")])
rownames(fdat) <- rownames(htseq_counts1)
esetRaw <- to.eSet(mat=htseq_counts1, pdat=sample_info, fdat=fdat)

## remove 'undetermined' samples
esetRaw <- esetRaw[,!is.na(pData(esetRaw)$Group)]

## FPKM (for limma)

## will use cufflinks fpkm for limma analysis
fpkm <- read.table(file=file.path(OMPATH,"data/zebrafish_cufflinks_counts_fpkm.txt"),
                   sep="\t", header=TRUE)[-1,] # remove first row (non-expression info)

## remove duplicates by taking those w/ highest variation
MAD <- apply(apply(fpkm[,-(1:2)],2,as.numeric),1,mad,na.rm=TRUE)
ord <- order(MAD,decreasing=TRUE)
fpkm <- fpkm[ord,]
fpkm1 <- fpkm[match(unique(fpkm[,"ENS_ID"]),fpkm[,"ENS_ID"]),]

rownames(fpkm1) <- fpkm1[,"ENS_ID"]   # use Ensemble IDs as rownames
fpkm1 <- fpkm1[,-(1:2)]               # .. remove column from data
fpkm1[,] <- apply(fpkm1,2,as.numeric) # .. and convert to numeric matrix

## start by creating an ExpressionSet from FPKM
cmn <- intersect( rownames(esetRaw), rownames(fpkm1) )
esetRaw <- esetRaw[cmn,]
esetFPKM.exprs <- fpkm1[match.nona(cmn, rownames(fpkm1)),
                        match.nona(colnames(esetRaw), colnames(fpkm1))]

all(colnames(esetFPKM.exprs)==colnames(esetRaw))
all(rownames(esetFPKM.exprs)==rownames(esetRaw))

esetFPKM <- to.eSet(mat=esetFPKM.exprs, pdat=pData(esetRaw), fdat=fData(esetRaw))
exprs(esetFPKM) <- log2(exprs(esetFPKM)+1)

## check one more time
all(featureNames(esetRaw)==featureNames(esetFPKM))
all(sampleNames(esetRaw)==sampleNames(esetFPKM))

## save both RAW and FPKM ExpressionSets
saveRDS(esetRaw, file=file.path(OMPATH,"data/zebrafish_htseq_raw_counts_eSet.RDS"))
saveRDS(esetFPKM,file=file.path(OMPATH,"data/zebrafish_cufflinks_counts_fpkm_eSet.RDS"))


## esetRaw <- readRDS(file=file.path(OMPATH,"data/zebrafish_htseq_raw_counts_eSet.RDS"))
