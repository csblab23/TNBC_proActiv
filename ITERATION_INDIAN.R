## RUN ITERATIONS ON INDIAN DATASET:
## NOTE: These Indian samples were processed using the chinese gtf file
## here we are taking REL PA > 0.05, 0.1, 0.15

setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf")

library("survival")
library("survminer")
library("proActiv")
library(dplyr)
library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
library(sva)
library(survival)

# load the metafile  
metadata <- as.data.frame(read_csv("./Indian_meta_TNBC.csv"))
metadata[1:5,1:5]
metadata=metadata[,-1]
dim(metadata) # 108  10
colnames(metadata)
table(metadata$Batch)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 15 
# 10  5  1 15 20 17  7  3  5  9  4  2  4  6 
table(metadata$Tumor)
# Tumor 
# 108
table(metadata$TNBC)
#  TNBC 
#  108 

#taking only TNBC sample from RG
subset_metadata = metadata %>% 
  filter(hospital == 'Rajiv Gandhi')
dim(subset_metadata) # 99  10
subset_metadata$`DFS_(DISEASE_FREE_SURVIVAL)` <- gsub("Days|days","",subset_metadata$`DFS_(DISEASE_FREE_SURVIVAL)`)
subset_metadata <- subset_metadata[!subset_metadata$`DFS_(DISEASE_FREE_SURVIVAL)`=="TREATMENT ONGOING",]
subset_metadata$`DFS_(DISEASE_FREE_SURVIVAL)`<- as.numeric(subset_metadata$`DFS_(DISEASE_FREE_SURVIVAL)`)
subset_metadata<- subset_metadata[subset_metadata$`DFS_(DISEASE_FREE_SURVIVAL)`>0,]
subset_metadata<- na.omit(subset_metadata)
dim(subset_metadata) # 91  10

# rename this column properly
dfs_col <- grep("DFS",colnames(subset_metadata))
colnames(subset_metadata)[dfs_col] <- "DFS_in_days"
dim(subset_metadata)  # 91  10

# reading RDS file
result2_indian <- readRDS('./result2_108TNBCs.rds')
dim(result2_indian) # 53910   108
rowData(result2_indian)
result2_promoter = rowData(result2_indian)
head(result2_promoter)
result2_promoter = as.data.frame(result2_promoter)
dim(result2_promoter) # 53910   8

# Read the feature counts file
count <- as.data.frame(read_csv("Indian_featurecounts.csv"))
dim(count) # 62700   110

# taking >=2 promoters gene
filtered_df <- result2_promoter %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep=''))
dim(filtered_df) #  24120     9

filtered_df <- filtered_df[, c('promoterId','geneId', 'ID')]
dim(filtered_df) #  24120   3
length(unique(filtered_df$geneId)) # 9315

table(table(filtered_df$geneId))
#    2    3    4    5    6    7    8    9   10   11   12   13   14   16   17   18   33 
# 6167 1941  697  274   95   59   36   13    8   10    7    1    1    3    1    1    1

# read the promoter counts
promoter_counts <- as.data.frame(result2_indian@assays@data@listData$promoterCounts)
length(intersect(filtered_df$promoterId, rownames(promoter_counts)))  # 24120
keep = intersect(filtered_df$promoterId, rownames(promoter_counts))
promoter_counts <- promoter_counts[keep,]
dim(promoter_counts) # 24120   108
promoter_counts[1:5,1:5]
colnames(promoter_counts) <- gsub("SJ.out", "", colnames(promoter_counts))
colnames(promoter_counts) <- gsub("^X", "", colnames(promoter_counts))
colnames(promoter_counts) <- gsub("_ME","",colnames(promoter_counts))
colnames(promoter_counts) <- gsub("\\.","-",colnames(promoter_counts))

subset_metadata$Sample <- gsub("_ME","",subset_metadata$Sample)

a=intersect(subset_metadata$Sample,colnames(promoter_counts))
length(a)  # 91
d = setdiff(subset_metadata$Sample,colnames(promoter_counts))
length(d)  #  0

promoter_counts <- promoter_counts[,colnames(promoter_counts) %in% subset_metadata$Sample]
promoter_counts <- promoter_counts[,subset_metadata$Sample]
all(colnames(promoter_counts)==subset_metadata$Sample)  # TRUE
dim(promoter_counts)  # 24120    91

# merge the DFs
p_counts <- merge( filtered_df, promoter_counts, by.x = "promoterId", by.y = 'row.names')
dim(p_counts) # 24120   94
p_counts[1:5,1:5]

# Select genes with atleast 2 promoters having count >5 or >10 (depending upon the iteration being run)
gene_lists <- list()
for (sample_col in colnames(p_counts)[4:ncol(p_counts)]) {
  cat("Processing sample:", sample_col, "\n")
  filtered_genes <- c()
  for (gene in unique(p_counts$geneId)) {
    gene_counts <- p_counts[p_counts$geneId == gene, sample_col]
    num_ids_gt_threshold <- sum(gene_counts > 5)
    if (num_ids_gt_threshold >= 2) {
      filtered_genes <- c(filtered_genes, gene)
    }
  }
  gene_lists[[sample_col]] <- filtered_genes
}

# Now take the genes that are present in atleast 60% of samples
n_samples <- length(gene_lists)
ensg_freq <- table(unlist(gene_lists))
threshold <- round(0.6 * n_samples)
ensg_60percent <- names(ensg_freq[ensg_freq >= threshold])
length(ensg_60percent) # 257
dim(filtered_df) #  24120     3
length(unique(filtered_df$geneId)) # 9315
filtered_df <- filtered_df[filtered_df$geneId %in% ensg_60percent,]
dim(filtered_df)  #  882   3
length(unique(filtered_df$geneId)) # 257
table(table(filtered_df$geneId))
#  2  3  4  5  6  7  8 11 12 18 33 
# 93 88 40 17  4  6  3  2  2  1  1 

# read the relative counts df now
Rl_counts <- as.data.frame(result2_indian@assays@data@listData$relativePromoterActivity)
dim(Rl_counts) # 53910   108
Rl_counts <- Rl_counts[complete.cases(Rl_counts),]
dim(Rl_counts) # 754 108
length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  # 418
keep = intersect(filtered_df$promoterId, rownames(Rl_counts))
Rl_counts <- Rl_counts[keep,]
dim(Rl_counts) # 418  108
Rl_counts[1:5,1:5]
colnames(Rl_counts) <- gsub("SJ.out", "", colnames(Rl_counts))
colnames(Rl_counts) <- gsub("^X","",colnames(Rl_counts))
colnames(Rl_counts) <- gsub("_ME","",colnames(Rl_counts))
colnames(Rl_counts) <- gsub("\\.","-",colnames(Rl_counts))

a=intersect(subset_metadata$Sample,colnames(Rl_counts))
length(a)  # 91
d=setdiff(subset_metadata$Sample,colnames(Rl_counts))
length(d)  # 0

Rl_counts <- Rl_counts[,colnames(Rl_counts) %in% a]
Rl_counts <- Rl_counts[,subset_metadata$Sample]
all(subset_metadata$Sample==colnames(Rl_counts))  # TRUE

## now apply loop: 
Rl_count_level <- c(0.05, 0.1, 0.15)
length(Rl_count_level) # 3
for (j in 1:length(Rl_count_level)){
  tobetaken_rel <- rowSums(Rl_counts > Rl_count_level[j]) >= dim(Rl_counts)[2]*0.1
  table(tobetaken_rel)
  Rl_counts=Rl_counts[tobetaken_rel,]
  dim(Rl_counts) #763 121
  length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  ## 763
  filtered_df <- filtered_df[filtered_df$promoterId %in% rownames(Rl_counts),]
  dim(filtered_df) #763   3
  length(unique(filtered_df$geneId)) #695
  table(table(filtered_df$geneId))
  ensg_freq_2 <- table(filtered_df$geneId)
  ensg_final <- names(ensg_freq_2[ensg_freq_2 >= 2])
  length(ensg_final) #88
  filtered_df <- filtered_df[filtered_df$geneId %in% ensg_final,]
  dim(filtered_df) #176   3
  length(unique(filtered_df$geneId)) #88
  table(table(filtered_df$geneId)) 
  
  count <- read_csv("Indian_featurecounts.csv") %>%
    as.data.frame() %>%
    column_to_rownames("Geneid") %>% 
    subset(., select = -1) %>% 
    filter(rownames(.) %in% filtered_df$geneId)
  dim(count) #88 420
  
  # taking same samples as metadata
  length(intersect(colnames(count), subset_metadata$Sample))
  keep = intersect(colnames(count), subset_metadata$Sample)
  count = count[, keep]
  count = count[,subset_metadata$Sample]
  all(subset_metadata$Sample==colnames(count))
  dim(count) # 62700    91
  cpm <- cpm(count)
  is.exprs <- rowSums(cpm>1) >= dim(count)[2]*0.1
  count <- count[is.exprs, ]
  x <- DGEList(counts=count )
  x <- calcNormFactors(x,method = "TMM")
  v <- voom(x, plot=F)
  vMat <- v$E
  
  # adjust batch effect:
  dim(subset_metadata) # 91  10
  table(subset_metadata$Batch)
  #  1  2  3  4  5  6  7  8  9 10 11 12 13 15 
  # 10  5  1 15 15 17  7  3  3  9  3  2  3  6 
  common_samples <- intersect(colnames(vMat), subset_metadata$Sample)
  vMat <- vMat[,common_samples]
  vMat <- vMat[,subset_metadata$Sample]
  all(colnames(vMat)==subset_metadata$Sample)
  
  subset_metadata$Batch = as.factor(subset_metadata$Batch)
  vMat_2 <- sva::ComBat(vMat, batch=subset_metadata$Batch)
  batch = subset_metadata$Batch
  dim(vMat_2) # 88 121
  
  length(intersect(filtered_df$geneId , rownames(vMat_2)))
  dim(filtered_df) #176   3
  p_info <- filtered_df[filtered_df$geneId %in% row.names(vMat_2),]
  dim(p_info) #176   3
  length(unique(p_info$geneId)) #88
  
  # read the absolute counts matrix
  Ab_counts <- as.data.frame(result2_indian@assays@data@listData$absolutePromoterActivity)
  Ab_counts <- Ab_counts[rownames(Ab_counts) %in% p_info$promoterId,]
  colnames(Ab_counts) <- gsub("SJ.out", "", colnames(Ab_counts))
  colnames(Ab_counts) <- gsub("^X","",colnames(Ab_counts))
  colnames(Ab_counts) <- gsub("_ME","",colnames(Ab_counts))
  colnames(Ab_counts) <- gsub("\\.","-",colnames(Ab_counts))
  
  a=intersect(subset_metadata$Sample,colnames(Ab_counts))
  length(a)  # 91
  d=setdiff(subset_metadata$Sample,colnames(Ab_counts))
  length(d)  # 0
  Ab_counts <- Ab_counts[,colnames(Ab_counts) %in% a,]
  Ab_counts <- Ab_counts[,subset_metadata$Sample]
  all(subset_metadata$Sample==colnames(Ab_counts))  # TRUE
  
  dim(Ab_counts) # 53910    91
  dim(p_info)
  length(intersect(p_info$promoterId , rownames(Ab_counts))) # 55873
  p_info <- p_info[ p_info$promoterId %in% rownames(Ab_counts) ,]
  dim(p_info) #176    3
  table(table(p_info$geneId))
  
  ensg_freq_3 <- table(p_info$geneId)
  ensg_final_1 <- names(ensg_freq_3[ensg_freq_3 >= 2])
  length(ensg_final_1) # 88
  p_info <- p_info[p_info$geneId %in% ensg_final_1,]
  dim(p_info) #176     3
  length(unique(p_info$geneId)) #88
  table(table(p_info$geneId)) 
  
  Ab_counts <- Ab_counts[rownames(Ab_counts) %in% p_info$promoterId,]
  Ab_counts <- merge( p_info, Ab_counts, by.x= 'promoterId', by.y = 'row.names')
  rownames(Ab_counts) <- Ab_counts$ID
  Ab_counts <- Ab_counts[, c(-1, -2, -3)]
  dim(Ab_counts) # 176  121
  
  str(subset_metadata)
  dim(subset_metadata) # 116   7
  subset_metadata_1 <- subset_metadata
  dim(subset_metadata)
  
  Ab_counts <- Ab_counts[, subset_metadata_1$Sample]
  all(subset_metadata_1$Sample==colnames(Ab_counts))
  dim(Ab_counts)# 176   113
  vMat_2 <- vMat_2[, subset_metadata_1$Sample]
  all(subset_metadata_1$Sample==colnames(vMat_2))
  dim(vMat_2) #88   113
  
  tobetaken <- rowSums(Ab_counts > 0) >= dim(subset_metadata_1)[1]*0.1
  Ab_counts=Ab_counts[tobetaken,]
  dim(Ab_counts) #176 113
  
  subset_metadata_1$`RECC/PROG_IF_ANY,` = as.factor(subset_metadata_1$`RECC/PROG_IF_ANY,`)
  subset_metadata_1$dfs_status <- ifelse(subset_metadata_1$`RECC/PROG_IF_ANY` == 'YES', 1, 0)
  
  
  
  surv_level = c(0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5)  
  for (l in 1:length(surv_level)){
    #function for performing matrix processing
    rownames(subset_metadata_1) = subset_metadata_1$Sample
    perform_analysis_for_matrix <- function(metadata, matrix) {
      mm <- merge(metadata, matrix, by = 0)
      rownames(mm) <- mm$Row.names
      mm <- mm[, -1]
    }
    prom_matrix <- perform_analysis_for_matrix(subset_metadata_1, t(Ab_counts))
    dim(prom_matrix) #113 185
    
    rownames(subset_metadata_1) = subset_metadata_1$Sample
    perform_analysis_for_matrix <- function(metadata, matrix) {
      mm <- merge(metadata, matrix, by = 0)
      rownames(mm) <- mm$Row.names
      mm <- mm[, -1]
    }
    
    gene_matrix <- perform_analysis_for_matrix(subset_metadata_1, t(vMat_2))
    dim(gene_matrix) #113  97
    
    prom_p_values <- list()
    gene_p_values <- list()
    
    #function for performing surv analysis
    perform_survival_analysis <- function(mm_check, column_name) {
      test_surv_df <- data.frame(Sample = mm_check$Sample, Event = mm_check$dfs_status, Time = mm_check$DFS_in_days)
      mm_check_subset <- mm_check[, c('Sample', column_name)]
      mm_check_subset <- mm_check_subset %>% arrange(desc(!!sym(column_name)))
      mm_check_subset <- mm_check_subset %>%
        mutate(Binary = ifelse(row_number() <= round(dim(mm_check_subset)[1]*(surv_level[l])), "high", "low"))
      merged_df <- merge(test_surv_df, mm_check_subset, by = "Sample")
      merged_df$Binary = as.factor(merged_df$Binary)
      fit <- coxph(Surv(merged_df$Time, merged_df$Event) ~ merged_df[,'Binary'], data = merged_df)
      temp_p <- summary(fit)$coefficients[, "Pr(>|z|)"]
      return(temp_p)
    }
    
    for (sample_col in colnames(prom_matrix)[13:ncol(prom_matrix)]) {
      prom_p_values[[sample_col]] <- perform_survival_analysis(prom_matrix, sample_col)
      cat("Promoter P-value for", sample_col, ":", prom_p_values[[sample_col]], "\n")
    }
    
    for (sample_col in colnames(gene_matrix)[13:ncol(gene_matrix)]) {
      gene_p_values[[sample_col]] <- perform_survival_analysis(gene_matrix, sample_col)
      cat("Gene P-value for", sample_col, ":", gene_p_values[[sample_col]], "\n")
    }
    
    fdr_p_values = p.adjust(unlist(prom_p_values), method = 'fdr')  
    prom_df <- data.frame(sample_ID = names(prom_p_values), p_value = unlist(prom_p_values), fdr_value = fdr_p_values)
    fdr_g_p_values = p.adjust(unlist(gene_p_values), method = 'fdr')  
    gene_df <- data.frame(sample_ID = names(gene_p_values), p_value = unlist(gene_p_values), fdr_values = fdr_g_p_values)
    
    prom_sig_result <- filter(prom_df, p_value < 0.05)
    prom_sig_result_adj <- filter(prom_df, fdr_value < 0.05)
    gene_sig_result <- filter(gene_df, p_value < 0.1)
    gene_sig_result_adj <- filter(gene_df, fdr_values < 0.1)
    gene_nonsig_result <- filter(gene_df, p_value > 0.1)
    gene_nonsig_result_adj <- filter(gene_df, fdr_values > 0.1)
    
    #intersection lists
    prom_sig_result_2 = prom_sig_result %>% separate(sample_ID,into =c("Promoter", "Gene_names2"),sep = "_")
    com_p_g_sig_pv = (intersect(prom_sig_result_2$Gene_names2,gene_sig_result$sample_ID))  #both_sig_p_value
    com_p_sig_g_nonsig_pv = (intersect(prom_sig_result_2$Gene_names2, gene_nonsig_result$sample_ID)) #by pv
    com_p_sig_pv_g_nonsig_fdr = (intersect(prom_sig_result_2$Gene_names2, gene_nonsig_result_adj$sample_ID))
    
    
    
    output <- data.frame(total_prom = dim(prom_df)[1], total_gene = dim(gene_df)[1],
                         sig_pr_pv = dim(prom_sig_result)[1], sig_pr_fdr = dim(prom_sig_result_adj)[1],
                         sig_g_pv = dim(gene_sig_result)[1], sig_g_fdr = dim(gene_sig_result_adj)[1],
                         nonsig_g_pv = dim(gene_nonsig_result)[1], nonsig_g_fdr = dim(gene_nonsig_result_adj),
                         both_sig_pv = length(com_p_g_sig_pv), p_sig_g_nonsig_pv = length(com_p_sig_g_nonsig_pv),
                         p_sig_g_nonsig_fdr = length(com_p_sig_pv_g_nonsig_fdr))
    
    if (j == 1){
      if (l == 1){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.05/surv_top_10%")
        cat("You are in 5_count/Rl_PA_0.1/surv_top_10% working directory")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      if (l ==2){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.05/surv_top_20%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      if (l == 3){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.05/surv_top_25%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
        
      }
      if ( l ==4){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.05/surv_top_30%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      if (l ==5){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.05/surv_top_35%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      if ( l == 6){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.05/surv_top_40%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      if (l == 7){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.05/surv_top_50%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      
    }
    if (j ==2){
      if (l == 1){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.1/surv_top_10%")
        cat("You are in 5_count/Rl_PA_0.01/surv_top_10% working directory")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l ==2){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.1/surv_top_20%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l == 3){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.1/surv_top_25%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
        
      }
      else if ( l ==4){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.1/surv_top_30%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l ==5){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.1/surv_top_35%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if ( l == 6){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.1/surv_top_40%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l == 7){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.1/surv_top_50%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      
    }
    if (j == 3){
      if (l == 1){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.15/surv_top_10%")
        cat("You are in 5_count/Rl_PA_0.15/surv_top_10% working directory")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l ==2){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.15/surv_top_20%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l == 3){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.15/surv_top_25%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
        
      }
      else if ( l ==4){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.15/surv_top_30%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l ==5){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.15/surv_top_35%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if ( l == 6){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.15/surv_top_40%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      else if (l == 7){
        setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf/10_counts/Rl_PA_0.15/surv_top_50%/")
        write.csv(gene_df, 'all_gene.csv')
        write.csv(prom_df, 'all_prom.csv')
        write.csv(prom_sig_result, "prom_sig_pvalue.csv")
        write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
        write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
        write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
        write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
        write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
        write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
        write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
        write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
        write.csv(output, 'output.csv')
      }
      
    }
    # if (j == 4){
    #   if (l == 1){
    #     setwd("C:/Users/Simran/Desktop/Promoter activity project/loop_pipeline_updated_March21/Iterations_may30/Indian_dataset/5_counts/Rl_PA_0.5/surv_top_10%")
    #     cat("You are in 5_counts/Rl_PA_0.5/surv_top_10% working directory")
    #     write.csv(gene_df, 'all_gene.csv')
    #     write.csv(prom_df, 'all_prom.csv')
    #     write.csv(prom_sig_result, "prom_sig_pvalue.csv")
    #     write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
    #     write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
    #     write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
    #     write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
    #     write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
    #     write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
    #     write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
    #     write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
    #     write.csv(output, 'output.csv')
    #   }
    #   if (l ==2){
    #     setwd("C:/Users/Simran/Desktop/Promoter activity project/loop_pipeline_updated_March21/Iterations_may30/Indian_dataset/5_counts/Rl_PA_0.5/surv_top_20%/")
    #     write.csv(gene_df, 'all_gene.csv')
    #     write.csv(prom_df, 'all_prom.csv')
    #     write.csv(prom_sig_result, "prom_sig_pvalue.csv")
    #     write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
    #     write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
    #     write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
    #     write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
    #     write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
    #     write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
    #     write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
    #     write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
    #     write.csv(output, 'output.csv')
    #   }
    #   if (l == 3){
    #     setwd("C:/Users/Simran/Desktop/Promoter activity project/loop_pipeline_updated_March21/Iterations_may30/Indian_dataset/5_counts/Rl_PA_0.5/surv_top_25%/")
    #     write.csv(gene_df, 'all_gene.csv')
    #     write.csv(prom_df, 'all_prom.csv')
    #     write.csv(prom_sig_result, "prom_sig_pvalue.csv")
    #     write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
    #     write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
    #     write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
    #     write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
    #     write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
    #     write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
    #     write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
    #     write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
    #     write.csv(output, 'output.csv')
    #     
    #   }
    #   if ( l ==4){
    #     setwd("C:/Users/Simran/Desktop/Promoter activity project/loop_pipeline_updated_March21/Iterations_may30/Indian_dataset/5_counts/Rl_PA_0.5/surv_top_30%/")
    #     write.csv(gene_df, 'all_gene.csv')
    #     write.csv(prom_df, 'all_prom.csv')
    #     write.csv(prom_sig_result, "prom_sig_pvalue.csv")
    #     write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
    #     write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
    #     write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
    #     write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
    #     write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
    #     write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
    #     write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
    #     write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
    #     write.csv(output, 'output.csv')
    #   }
    #   if (l ==5){
    #     setwd("C:/Users/Simran/Desktop/Promoter activity project/loop_pipeline_updated_March21/Iterations_may30/Indian_dataset/5_counts/Rl_PA_0.5/surv_top_35%/")
    #     write.csv(gene_df, 'all_gene.csv')
    #     write.csv(prom_df, 'all_prom.csv')
    #     write.csv(prom_sig_result, "prom_sig_pvalue.csv")
    #     write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
    #     write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
    #     write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
    #     write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
    #     write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
    #     write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
    #     write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
    #     write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
    #     write.csv(output, 'output.csv')
    #   }
    #   if ( l == 6){
    #     setwd("C:/Users/Simran/Desktop/Promoter activity project/loop_pipeline_updated_March21/Iterations_may30/Indian_dataset/5_counts/Rl_PA_0.5/surv_top_40%/")
    #     write.csv(gene_df, 'all_gene.csv')
    #     write.csv(prom_df, 'all_prom.csv')
    #     write.csv(prom_sig_result, "prom_sig_pvalue.csv")
    #     write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
    #     write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
    #     write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
    #     write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
    #     write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
    #     write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
    #     write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
    #     write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
    #     write.csv(output, 'output.csv')
    #   }
    #   if (l == 7){
    #     setwd("C:/Users/Simran/Desktop/Promoter activity project/loop_pipeline_updated_March21/Iterations_may30/Indian_dataset/5_counts/Rl_PA_0.5/surv_top_50%/")
    #     write.csv(gene_df, 'all_gene.csv')
    #     write.csv(prom_df, 'all_prom.csv')
    #     write.csv(prom_sig_result, "prom_sig_pvalue.csv")
    #     write.csv(prom_sig_result_adj, "prom_sig_fdr.csv")
    #     write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
    #     write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
    #     write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
    #     write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')
    #     write.csv(com_p_g_sig_pv, 'com_p_g_sig_pv.csv')
    #     write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')
    #     write.csv(com_p_sig_pv_g_nonsig_fdr, 'com_p_sig_pv_g_nonsig_fdr.csv')
    #     write.csv(output, 'output.csv')
    #   }
      
    }
    
  setwd("~/Promoter_activity/proActiv/indian_dataset/iteration_28may/with_new_gtf")
    
  }
  




####################################










