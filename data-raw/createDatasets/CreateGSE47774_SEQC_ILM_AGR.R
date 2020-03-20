## DATA DOWNLOADED FROM
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47774
##
## Expression Matrix File (in compressed txt format):
## https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE47774&format=file&file=GSE47774%5FSEQC%5FILM%5FAGR%2Etxt%2Egz
##
## Data description from: http://www.nature.com/nbt/journal/v32/n9/full/nbt.2931.html
##
## SEQC data set. The third phase of the MicroArray Quality Control
## (MAQC) project, also known as the Sequencing Quality Control (SEQC)
## project, aims to assess the technical performance of
## high-throughput sequencing platforms by generating benchmarking
## data sets. The design includes four different sample types, namely
## samples A, B, C and D. Sample A is Stratagene's universal human
## reference (UHR) RNA; sample B is Ambion's human brain reference
## RNA; samples C and D are mixes of samples A and B, in a 3:1 and 1:3
## ratio, respectively. The four reference samples were sent to
## several sequencing centers around the world and sequenced using
## different platforms. Here, we focus on sample A and sample B
## sequenced at the Australian Genome Research Facility (AGRF) using
## the Illumina HiSeq2000. Four libraries were prepared for each of
## sample A and B and multiplex pools of the resulting 8 barcoded
## libraries were sequenced in 8 lanes of 2 flow-cells, yielding a
## total of 16 (technical) replicates per library and 64 replicates
## per sample type. Prior to library preparation, Ambion ERCC ExFold
## RNA Spike-in Control Mix 1 and Mix 2 were added to sample A and
## sample B RNA, respectively, in a proportion of 50 μl per 2,500 μl
## of total RNA. The data consist of an average of 10 million 100-bp
## paired-end reads per sample.

require(Biobase)
OMPATH <- Sys.getenv("OMPATH")

## read the txt file (after uncompressing it)
system( paste("gunzip ",OMPATH,"/data/GSE47774_SEQC_ILM_AGR.txt.gz",sep="") )
SEQC <- read.delim(file.path(OMPATH,"data/GSE47774_SEQC_ILM_AGR.txt"),row.names=1)
system( paste("gzip ",OMPATH,"/data/GSE47774_SEQC_ILM_AGR.txt",sep="") )

## create the sample annotation table
##
ANN <- as.data.frame(t(sapply(gsub("SEQC_ILM_AGR_","",colnames(SEQC)),
                              function(Z) unlist(strsplit(Z,"_")))))
rownames(ANN) <- colnames(SEQC)
colnames(ANN) <- c("type","library","lane","barcode","flowcell")

## retrieve probes' annotation (removed 'ensemble_id', 'chromosome',
## etc. from attributes to minimize duplicates)

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
GS <- getBM(filters="refseq_mrna",
            attributes= c("refseq_mrna","hgnc_symbol","description"),
            values=rownames(SEQC),
            mart= mart)

## quick inspection of the results
sum(tmp <- duplicated(GS[,1]))         # 30
GS[grep(GS[which(tmp)[1],1],GS[,1]),]  # NM_000344 --> {SMN2, SMN1}

## create the fData, filling in the non-empty entries returned by getBM
##
fdata <- data.frame(refseq=rownames(SEQC),hgnc_symbol="",description="",stringsAsFactors=FALSE,row.names=1)
matchIDX <- match(rownames(SEQC),GS[,1])
fdata[!is.na(matchIDX),] <- GS[matchIDX[!is.na(matchIDX)],c("hgnc_symbol","description")]
fdata[grep("ERCC",rownames(fdata)),"description"] <- "spike-in"

## finally, create the ExpressionSet
ESET <- ExpressionSet(assayData=as.matrix(SEQC),
                      phenoData=AnnotatedDataFrame(ANN),
                      featureData=AnnotatedDataFrame(fdata))
saveRDS(ESET,file=file.path(OMPATH,"data/GSE47774_SEQC_ILM_AGR.rawdata.eSet.RDS"))

## Also, create a version restricted to entries w/ non-empty gene symbols,
## and type A and B samples only (drop empty factor levels from the pData)

ESET1 <- ESET[fData(ESET)[,"hgnc_symbol"]!="", ESET$type %in% c("A","B")]
for ( i in 1:ncol(pData(ESET1)) ) { pData(ESET1)[,i] <- droplevels(pData(ESET1)[,i]) }
saveRDS(ESET1,file=file.path(OMPATH,"data/GSE47774_SEQC_ILM_AGR.ABonly.wGeneSymbol.eSet.RDS"))
