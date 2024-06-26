#AIM: to split the datasets into training and testing dataset = CARET package 

setwd("~/proActiv/univariate_surv_proActiv")

# Load required packages
library("survival")
library("survminer")
library("proActiv")
library(data.table)
library(dplyr)
library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(sva)

#loading metafile
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run

tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor')
dim(tumor_metadata) #360  99

colnames(tumor_metadata)
surv_metadata =tumor_metadata[tumor_metadata$RFS_time_Months > 0,]
dim(surv_metadata) #358  99
surv_metadata=surv_metadata[,c("Project_ID", "Run", "batch", "Tumor_Normal", "Tumor_Normal1", "RFS_Status", "RFS_time_Days", "RFS_time_Months", 'CNA_Subtype', 'Grade')]
dim(surv_metadata) #358  10

filter_data <- surv_metadata %>%
  mutate(CNA_Subtype = replace_na(CNA_Subtype, 'unknown'),
         Grade = replace_na(Grade, 'unknown'))

dim(filter_data) #358  10
table(filter_data$CNA_Subtype)
# Chr12p13_amp Chr13q34_amp Chr20q13_amp  Chr8p21_del  Chr9p23_amp      Low_CIN      unknown 
# 46           25           30           48           41          110           58 
table(filter_data$Grade)
# 2  2 to 3       3 unknown 
# 64      29     231      34 

str(filter_data)
filter_data$CNA_Subtype = as.factor(filter_data$CNA_Subtype)
filter_data$Grade = as.factor(filter_data$Grade)
str(filter_data)


#dividing the data in train-test split using caret package
library(caret)

filter_data_2 <- filter_data %>%  mutate(combo = paste(filter_data$CNA_Subtype, filter_data$Grade, sep = '___') )

str(filter_data_2)
filter_data_2$combo <- as.factor(filter_data_2$combo)
str(filter_data_2)
levels(filter_data_2$combo) #27 levels
table(filter_data_2$combo)

#to check how many levels, in filter_data_2 have 


set.seed(123)
inTrain <- createDataPartition(
  y = filter_data_2$combo,
  p = .70, ## Proportion for the training set
  list = FALSE
)

str(inTrain)
training <- filter_data_2[ inTrain,]
testing  <- filter_data_2[-inTrain,]
dim(training) #263  11 #split may vary slightly due to the random sampling process --- number of levels that have only 1 sample would vary the split number 
dim(testing) #95 11

#validation
table(training$CNA_Subtype)
table(testing$CNA_Subtype)
table(training$Grade)
table(testing$Grade)


prop.table(table(training$CNA_Subtype))
# Chr12p13_amp Chr13q34_amp Chr20q13_amp  Chr8p21_del  Chr9p23_amp      Low_CIN      unknown 
# 0.13307985   0.07224335   0.08745247   0.13307985   0.11406844   0.30038023   0.15969582 
prop.table(table(testing$CNA_Subtype))
# Chr12p13_amp Chr13q34_amp Chr20q13_amp  Chr8p21_del  Chr9p23_amp      Low_CIN      unknown 
# 0.11578947   0.06315789   0.07368421   0.13684211   0.11578947   0.32631579   0.16842105 


prop.table(table(training$Grade))
# 2     2 to 3          3    unknown 
# 0.18250951 0.08365019 0.62357414 0.11026616 
prop.table(table(testing$Grade))
# 2     2 to 3          3    unknown 
# 0.16842105 0.07368421 0.70526316 0.05263158 


#now taking the training dataset for calculating cut-off, 
dim(training) #263  11
training[1:10, 1:11]



write.csv(x = training, file = 'training_dataset_v2.csv')
write.csv(x = testing, file = 'testing_dataset_v2.csv')
