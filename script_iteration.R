#for running the iterations

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

# Read input files
result2 <- readRDS('result2_v2_updated.rds')
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
count <- as.data.frame(read_csv("Featurecounts_combined_fuscc.csv"))

# Data preprocessing
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run

tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor')

#taking >=2 promoters gene
result2_promoter = rowData(result2) 
filtered_df <- result2_promoter %>%
  as.data.frame() %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

filtered_df <- filtered_df[, c('promoterId','geneId', 'ID')]
dim(filtered_df) #24120     3
length(unique(filtered_df$geneId)) #9315

table(table(filtered_df$geneId))
# 2    3    4    5    6    7    8    9   10   11   12   13   14   16   17   18   33 
# 6167 1941  697  274   95   59   36   13    8   10    7    1    1    3    1    1    1 

promoter_counts <- as.data.frame(result2@assays@data@listData$promoterCounts)
length(intersect(filtered_df$promoterId, rownames(promoter_counts)))  ## 24120
keep = intersect(filtered_df$promoterId, rownames(promoter_counts))
promoter_counts <- promoter_counts[keep,]
dim(promoter_counts) #24120   360
colnames(promoter_counts) <- gsub("_SJ.out", "", colnames(promoter_counts))
#merge the DFs
p_counts <- merge( filtered_df, promoter_counts, by.x = "promoterId", by.y = 'row.names')
dim(p_counts) ##24120   363

count_level = c(5, 10)
length(count_level)
for (i in 1:length(count_level)){
  gene_lists <- list()
  
  for (sample_col in colnames(p_counts)[4:ncol(p_counts)]) {
    cat("Processing sample:", sample_col, "\n")
    filtered_genes <- c()
    for (gene in unique(p_counts$geneId)) {
      gene_counts <- p_counts[p_counts$geneId == gene, sample_col]
      num_ids_gt_threshold <- sum(gene_counts > count_level[i])
      if (num_ids_gt_threshold >= 2) {
        filtered_genes <- c(filtered_genes, gene)
      }
    }
    gene_lists[[sample_col]] <- filtered_genes
  }
  # #for taking genes that are present in atleast 60% of samples ##---------------------------LEVEL-03
  n_samples <- length(gene_lists)
  ensg_freq <- table(unlist(gene_lists))
  threshold <- 0.6 * n_samples
  ensg_60percent <- names(ensg_freq[ensg_freq >= threshold])
  length(ensg_60percent) 
  dim(filtered_df) 
  length(unique(filtered_df$geneId)) 
  filtered_df <- filtered_df[filtered_df$geneId %in% ensg_60percent,]
  dim(filtered_df) 
  length(unique(filtered_df$geneId))
  table(table(filtered_df$geneId))
  
  Rl_counts <- as.data.frame(result2@assays@data@listData$relativePromoterActivity) ##----------LEVEL-04
  dim(Rl_counts) 
  Rl_counts <- Rl_counts[complete.cases(Rl_counts),]
  dim(Rl_counts) 
  length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  
  keep = intersect(filtered_df$promoterId, rownames(Rl_counts))
  Rl_counts <- Rl_counts[keep,]
  dim(Rl_counts) 
  colnames(Rl_counts) <- gsub("_SJ.out", "", colnames(Rl_counts))
  
  Rl_count_level <- c(0.2, 0.25, 0.3, 0.5)
  length(Rl_count_level)
  for (j in 1:length(Rl_count_level)){
    tobetaken_rel <- rowSums(Rl_counts > Rl_count_level[j]) >= dim(Rl_counts)[2]*0.1
    table(tobetaken_rel)
    Rl_counts=Rl_counts[tobetaken_rel,]
    dim(Rl_counts) 
    length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  
    filtered_df <- filtered_df[filtered_df$promoterId %in% rownames(Rl_counts),]
    dim(filtered_df) 
    length(unique(filtered_df$geneId)) 
    table(table(filtered_df$geneId))
    ensg_freq_2 <- table(filtered_df$geneId)
    ensg_final <- names(ensg_freq_2[ensg_freq_2 >= 2])
    length(ensg_final) #243
    filtered_df <- filtered_df[filtered_df$geneId %in% ensg_final,]
    dim(filtered_df) #486   3  #696   3
    length(unique(filtered_df$geneId)) #243 #348 
    table(table(filtered_df$geneId)) 
    
    
    #now we will pre-process feature count matrix here only, why here??
    #because we have to take genes that are present in the filtered step till here
    count <- read_csv("Featurecounts_combined_fuscc.csv") %>%
      as.data.frame() %>%
      column_to_rownames("Geneid") %>% 
      subset(., select = -1) %>% 
      filter(rownames(.) %in% filtered_df$geneId)
    cpm <- cpm(count)
    is.exprs <- rowSums(cpm>1) >= dim(tumor_metadata)[1]*0.1
    count <- count[is.exprs, ]
    x <- DGEList(counts=count )
    x <- calcNormFactors(x,method = "TMM")
    v <- voom(x, plot=F)
    vMat <- v$E
    vMat_2 <- sva::ComBat(vMat, batch=metadata$batch)
    vMat_2 <- vMat_2[, tumor_metadata$Run] #removing normal sample
    dim(vMat_2) 
    
    p_info <- filtered_df[filtered_df$geneId %in% row.names(vMat_2),]
    dim(p_info) 
    length(unique(p_info$geneId))
    
    #taking Abs PA
    Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
    Ab_counts <- Ab_counts[rownames(Ab_counts) %in% p_info$promoterId,]
    colnames(Ab_counts) <- gsub("_SJ.out", "", colnames(Ab_counts))
    Ab_counts <- merge( p_info, Ab_counts, by.x= 'promoterId', by.y = 'row.names')
    rownames(Ab_counts) <- Ab_counts$ID
    Ab_counts <- Ab_counts[, c(-1, -2, -3)]
    #TAKING SAMPLES THAT HAVE RFS_time_months > 0 
    surv_metadata =tumor_metadata[tumor_metadata$RFS_time_Months > 0,]
    surv_metadata=surv_metadata[,c("Project_ID", "Run", "batch", "Tumor_Normal", "Tumor_Normal1", "RFS_Status", "RFS_time_Days", "RFS_time_Months")]
    dim(surv_metadata) 
    Ab_counts <- Ab_counts[, surv_metadata$Run]
    vMat_2 <- vMat_2[, surv_metadata$Run]
    dim(Ab_counts) 
    tobetaken <- rowSums(Ab_counts > 0) >= dim(tumor_metadata)[1]*0.1
    Ab_counts=Ab_counts[tobetaken,]
    dim(Ab_counts) 
                                                           ####--------------------------------LEVEL- 07
    surv_level = c(0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5)  
    for (l in 1:length(surv_level)){
      #function for performing surv analysis
      perform_survival_analysis <- function(mm_check, column_name) {
        test_surv_df <- data.frame(Run = mm_check$Run, Event = mm_check$RFS_Status, Time = mm_check$RFS_time_Days)
        mm_check_subset <- mm_check[, c('Run', column_name)]
        mm_check_subset <- mm_check_subset %>% arrange(desc(!!sym(column_name)))
        mm_check_subset <- mm_check_subset %>%
          mutate(Binary = ifelse(row_number() <= round(dim(mm_check_subset)[1]*(surv_level[l])), "high", "low"))
        merged_df <- merge(test_surv_df, mm_check_subset, by = "Run")
        merged_df$Binary = as.factor(merged_df$Binary)
        fit <- coxph(Surv(merged_df$Time, merged_df$Event) ~ merged_df[,'Binary'], data = merged_df)
        temp_p <- summary(fit)$coefficients[, "Pr(>|z|)"]
        return(temp_p)
      }
      
      #function for performing matrix processing
      perform_analysis_for_matrix <- function(metadata, matrix) {
        mm <- merge(metadata, matrix, by = 0)
        rownames(mm) <- mm$Row.names
        mm <- mm[, -1]
      }
      
      prom_p_values <- list()
      gene_p_values <- list()
      
      prom_matrix <- perform_analysis_for_matrix(surv_metadata, t(Ab_counts))
      gene_matrix <- perform_analysis_for_matrix(surv_metadata, t(vMat_2))
      
      for (sample_col in colnames(prom_matrix)[9:ncol(prom_matrix)]) {
        prom_p_values[[sample_col]] <- perform_survival_analysis(prom_matrix, sample_col)
        cat("Promoter P-value for", sample_col, ":", prom_p_values[[sample_col]], "\n")
      }
      
      for (sample_col in colnames(gene_matrix)[9:ncol(gene_matrix)]) {
        gene_p_values[[sample_col]] <- perform_survival_analysis(gene_matrix, sample_col)
        cat("Gene P-value for", sample_col, ":", gene_p_values[[sample_col]], "\n")
      }
      
      fdr_p_values = p.adjust(unlist(prom_p_values), method = 'fdr')  
      prom_df <- data.frame(sample_ID = names(prom_p_values), p_value = unlist(prom_p_values), fdr_value = fdr_p_values)
      fdr_g_p_values = p.adjust(unlist(gene_p_values), method = 'fdr')  
      gene_df <- data.frame(sample_ID = names(gene_p_values), p_value = unlist(gene_p_values), fdr_values = fdr_g_p_values)
      
      #condition for printing out files
      if (i == 1){
        setwd("~/proActiv/iteration/5_counts")
        cat("You are in 5_count working directory")
        if (j == 1){
          setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2")
          cat("you are in 5_count/Rl_PA_0.2 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2/surv_top_10%")
            cat("You are in 5_count/Rl_PA_0.2/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.2/surv_top_50%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
         
        } 
        else if (j == 2){
          setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/")
          cat("you are in 5_count/Rl_PA_0.25 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/surv_top_10%")
            cat("You are in 5_count/Rl_PA_0.25/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.25/surv_top_50%/")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          
        }
        else if (j == 3){
          setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/")
          cat("you are in 5_count/Rl_PA_0.3 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/surv_top_10%")
            cat("You are in 5_count/Rl_PA_0.3/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.3/surv_top_50%/")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          
        }
        else {
          setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/")
          cat("you are in 5_count/Rl_PA_0.5 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/surv_top_10%")
            cat("You are in 5_count/Rl_PA_0.5/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/5_counts/Rl_PA_0.5/surv_top_50%/")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
        } 
          
      } 
      else {
        setwd("~/proActiv/iteration/10_counts")
        cat("You are in 10_count working directory")
        if (j == 1){
          setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2")
          cat("you are in 10_count/Rl_PA_0.2 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2/surv_top_10%")
            cat("You are in 10_count/Rl_PA_0.2/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.2/surv_top_50%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          
        } 
        else if (j == 2){
          setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/")
          cat("you are in 10_count/Rl_PA_0.25 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/surv_top_10%")
            cat("You are in 10_count/Rl_PA_0.25/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.25/surv_top_50%/")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          
        }
        else if (j == 3){
          setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/")
          cat("you are in 10_count/Rl_PA_0.3 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/surv_top_10%")
            cat("You are in 10_count/Rl_PA_0.3/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.3/surv_top_50%/")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          
        }
        else {
          setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/")
          cat("you are in 10_count/Rl_PA_0.5 working directory")
          if (l == 1){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/surv_top_10%")
            cat("You are in 10_count/Rl_PA_0.5/surv_top_10% working directory")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 2){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/surv_top_20%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 3){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/surv_top_25%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if ( l == 4){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/surv_top_30%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 5){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/surv_top_35%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 6){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/surv_top_40%")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
          else if (l == 7){
            setwd("~/proActiv/iteration/10_counts/Rl_PA_0.5/surv_top_50%/")
            write.csv(gene_df, 'all_gene.csv')
            write.csv(prom_df, 'all_prom.csv')
          }
        } 
        
      } 
    }
    setwd("~/proActiv/univariate_surv_proActiv")
  }
} 













