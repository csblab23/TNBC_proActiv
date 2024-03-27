#ProActiv pipeline

#input files: Metafile, SJ.out files, gtf files, 

setwd("~/Promoter_activity/proActiv/Pipeline_updated_March21")
library("survival")
library("survminer")
library("proActiv")
library(dplyr)
library(tidyverse)
library(readxl)

# load the metafile  
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata[1:5,1:5]
metadata=metadata[,-1]  # removing the extra s.no. column
dim(metadata) # 448  98
colnames(metadata)

metadata$Tumor_Normal1=metadata$Tumor_Normal
dim(metadata)  #448  99

metadata =metadata %>% 
  mutate(Tumor_Normal1 = case_when(Tumor_Normal=="TumorTissue"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.2018"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.rep"  ~ "Tumor",
                                   Tumor_Normal=="PairedNormalTissue"  ~ "AdjNormal",
                                   Tumor_Normal=="PairedNormalTissue.rep"  ~ "AdjNormal"))
table(metadata$Tumor_Normal1)
rownames(metadata) <- metadata$Run
metadata[1:5,1:5]

## Prepare the input data
# list the STAR junction files as input
files <- list.files(path = './all_SJ_files/', full.names = TRUE)  # 448 files
files2=as.data.frame(files)
head(files2,2)
files2$files=gsub("./all_SJ_files//","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
dim(files2) #448   1
files2$original = files
head(files2,4)

# Check if all samples are intersecting
length(intersect(files2$files, metadata$Run)) #448

# files2_new = files2[files2$files %in% metadata$Run,]
# dim(files2_new) #448   2

## Prepare promoter annotation
gtf.file <-('/home/dell/Promoter_activity/human_ref_genome/gencode.v44.basic.annotation.gtf')
#names(GenomeInfoDb::genomeStyles()) # To identify species argument to be used
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Before running ProActiv make sure to match the order of samples in metafile and files2 
metadata[1:5,1:5]
head(files2,2)
metadata = metadata[files2$files,]
all(rownames(metadata) %in% files2$files)  # TRUE
all(rownames(metadata)== files2$files)

allfilenames = files2$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE
class(allfilenames)
length(allfilenames) #448
head(allfilenames,4)
metadata[1:5,1:5]
tail(metadata)
head(allfilenames)
tail(allfilenames)

# Define the conditions
table(metadata$Tumor_Normal1)
# AdjNormal     Tumor 
#       88       360
class(metadata$Tumor_Normal1)
#metadata$Tumor_Normal1 <- as.factor(metadata$Tumor_Normal1)    # not converting to factor here
conditions <- metadata$Tumor_Normal1
promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = conditions)
show(result)
# class: SummarizedExperiment 
# dim: 77915 448 

rowData(result)
result_promoter <- rowData(result)
head(result_promoter)
result_promoter = as.data.frame(result_promoter)
dim(result_promoter) # 77915    14

#total promoters identify 
length(unique(result_promoter$promoterId)) # 77915

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
dim(result2) # 53910   448
rowData(result2)
result2_promoter = rowData(result2)
head(result2_promoter)
result2_promoter = as.data.frame(result2_promoter)
dim(result2_promoter) # 53910    14






