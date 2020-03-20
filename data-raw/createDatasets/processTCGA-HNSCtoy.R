## Load required libraries
library(Biobase)

## Load the processed RNASeq expression set object for all Head and Neck Squamous Cell Carcinoma (HNSC) TCGA samples

## Path to file on scc
src_dir <- "/restricted/projectnb/montilab-p/CBMrepositoryData/TCGA/ESets/"

## File corresponding to the latest processed ESet
filename <- "HNSC_2015_11_01_ES.rds"

HNSC.ES <- readRDS(file.path(src_dir,filename))

## Define some imoportant phenotypic data we want to retain
pheno_labels_to_keep <- c("bcr_patient_barcode","bcr_sample_barcode","tissue_type","my_grade","my_stage","patient.anatomic_neoplasm_subdivision")

## Take a look at the phenoData for these fields
head(pData(HNSC.ES)[,pheno_labels_to_keep])

## Reduce the phenotypic information, keeping only these fields of interest
phenoData(HNSC.ES) <- AnnotatedDataFrame(pData(HNSC.ES)[,pheno_labels_to_keep])

## Save to a new file

## Folder where output should be stored
out_dir <- "/restricted/projectnb/dynamicgep/multiportal/ESets/"
saveRDS(HNSC.ES, file.path(out_dir,"HNSC_RNASeq_ES.rds"))

## Down-sample from the eSet to make a smaller (toy) ESet
## Let's use 25 samples total
## Here, we'd like an even representation of tumor grades, so down-sample from each category
## Use a particular seed so that you get a certain 'random' each time
set.seed(123)
samples.to.keep <- c(sample(as.character(HNSC.ES$bcr_sample_barcode)[HNSC.ES$my_grade %in% "AN"],5),
                     sample(as.character(HNSC.ES$bcr_sample_barcode)[HNSC.ES$my_grade %in% "g1"],5),
                     sample(as.character(HNSC.ES$bcr_sample_barcode)[HNSC.ES$my_grade %in% "g2"],5),
                     sample(as.character(HNSC.ES$bcr_sample_barcode)[HNSC.ES$my_grade %in% "g3"],5),
                     sample(as.character(HNSC.ES$bcr_sample_barcode)[HNSC.ES$my_grade %in% "g4"],5))

toy.ES <- HNSC.ES[,samples.to.keep]

## One last filtering step to remove genes with 0 expression across all samples
toy.ES <- toy.ES[rowSums(exprs(toy.ES))!=0,]

## Save to file
saveRDS(toy.ES,file.path(out_dir,"HNSC_RNASeq_toy_ES.rds"))
