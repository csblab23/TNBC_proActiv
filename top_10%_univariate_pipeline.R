#univariate survival analysis: Max-stat, Union


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

# Read input files
result2 <- readRDS('result2.rds')
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
#around 6167 genes have around 2 promoters rest have more than 2 promoters
# 2    3    4    5    6    7    8    9   10   11   12   13   14   16   17   18   33 
# 6167 1941  697  274   95   59   36   13    8   10    7    1    1    3    1    1    1 

promoter_counts <- as.data.frame(result2@assays@data@listData$promoterCounts)
length(intersect(filtered_df$promoterId, rownames(promoter_counts)))  ## 24120
keep = intersect(filtered_df$promoterId, rownames(promoter_counts))
promoter_counts <- promoter_counts[keep,]
dim(promoter_counts) #24120   448
#taking tumor sample only
colnames(promoter_counts) <- gsub("_SJ.out", "", colnames(promoter_counts))
promoter_counts <- promoter_counts[,tumor_metadata$Run]
#merge the DFs
p_counts <- merge( filtered_df, promoter_counts, by.x = "promoterId", by.y = 'row.names')
dim(p_counts) ##24120   363

#selection criteria
gene_lists <- list()
for (sample_col in colnames(p_counts)[4:ncol(p_counts)]) {
  cat("Processing sample:", sample_col, "\n")
  filtered_genes <- c()
  for (gene in unique(p_counts$geneId)) {
    gene_counts <- p_counts[p_counts$geneId == gene, sample_col]
    num_ids_gt_threshold <- sum(gene_counts > 10)
    if (num_ids_gt_threshold >= 2) {
      filtered_genes <- c(filtered_genes, gene)
    }
  }
  gene_lists[[sample_col]] <- filtered_genes
}

#all_common_genes <- gene_lists[[1]] #store intersection of genes
# all_common_genes1 <-  gene_lists[[1]] #store union of genes

# #for taking genes that are present in atleast 60% of samples
n_samples <- length(gene_lists)
ensg_freq <- table(unlist(gene_lists))
threshold <- 0.6 * n_samples
ensg_60percent <- names(ensg_freq[ensg_freq >= threshold])
length(ensg_60percent) #828 #1278
dim(filtered_df) #24120     3
length(unique(filtered_df$geneId)) #9315
filtered_df <- filtered_df[filtered_df$geneId %in% ensg_60percent,]
dim(filtered_df) #2556    3   #3884    3
length(unique(filtered_df$geneId)) #828 #1278
table(table(filtered_df$geneId))
# 2   3   4   5   6   7   8   9  10  11  12  16  33 
# 390 247  94  40  20  14  11   3   2   3   2   1   1 

# for (i in 2:length(gene_lists)) {
#   #all_common_genes <- intersect(all_common_genes, gene_lists[[i]])
#   all_genes <- union(all_common_genes1, gene_lists[[i]])
# }

#takinng the genes that have promoter with relative activity > 0.5 in at least %10 of the tumor samples 
Rl_counts <- as.data.frame(result2@assays@data@listData$relativePromoterActivity)
dim(Rl_counts) #53910   448
Rl_counts <- Rl_counts[complete.cases(Rl_counts),]
dim(Rl_counts) #8355  448
length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  #2174 #3016
keep = intersect(filtered_df$promoterId, rownames(Rl_counts))
Rl_counts <- Rl_counts[keep,]
dim(Rl_counts) #2174  448
#taking tumor sample only
colnames(Rl_counts) <- gsub("_SJ.out", "", colnames(Rl_counts))
Rl_counts <- Rl_counts[,tumor_metadata$Run]
dim(Rl_counts) #2174  360

#TAKING only the promoters with relative activity > 0.5 in at least %10 of the tumor samples 
tobetaken_rel <- rowSums(Rl_counts > 0.5) >= dim(Rl_counts)[2]*0.1
table(tobetaken_rel)
# #tobetaken_rel
# FALSE  TRUE 
# 1393   781 
Rl_counts=Rl_counts[tobetaken_rel,]
dim(Rl_counts) #781 360
length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  ## 781
filtered_df <- filtered_df[filtered_df$promoterId %in% rownames(Rl_counts),]
dim(filtered_df) #781   3
length(unique(filtered_df$geneId)) #538 #791
table(table(filtered_df$geneId))
# 1   2 
# 295 243 
#have to take in consideration of genes that have more than 2 promoters after rel PA filter
ensg_freq_2 <- table(filtered_df$geneId)
ensg_final <- names(ensg_freq_2[ensg_freq_2 >= 2])
length(ensg_final) #243
filtered_df <- filtered_df[filtered_df$geneId %in% ensg_final,]
dim(filtered_df) #486   3  #696   3
length(unique(filtered_df$geneId)) #243 #348 
table(table(filtered_df$geneId)) 

#pre-processing on feature count matrix
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
#batch adjustment
library(sva)
vMat_2 <- sva::ComBat(vMat, batch=metadata$batch)
vMat_2 <- vMat_2[, tumor_metadata$Run] #removing normal sample
dim(vMat_2) #242 360

p_info <- filtered_df[filtered_df$geneId %in% row.names(vMat_2),]
dim(p_info) #484   3
length(unique(p_info$geneId)) #242

#taking Abs PA
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
Ab_counts <- Ab_counts[rownames(Ab_counts) %in% p_info$promoterId,]
colnames(Ab_counts) <- gsub("_SJ.out", "", colnames(Ab_counts))
Ab_counts <- Ab_counts[,tumor_metadata$Run]
Ab_counts <- merge( p_info, Ab_counts, by.x= 'promoterId', by.y = 'row.names')
rownames(Ab_counts) <- Ab_counts$ID
Ab_counts <- Ab_counts[, c(-1, -2, -3)]

#TAKING SAMPLES THAT HAVE RFS_time_months > 0 
surv_metadata =tumor_metadata[tumor_metadata$RFS_time_Months > 0,]
surv_metadata=surv_metadata[,c("Project_ID", "Run", "batch", "Tumor_Normal", "Tumor_Normal1", "RFS_Status", "RFS_time_Days", "RFS_time_Months")]
dim(surv_metadata) #358   8

Ab_counts <- Ab_counts[, surv_metadata$Run]
vMat_2 <- vMat_2[, surv_metadata$Run]
dim(Ab_counts) #484 358
tobetaken <- rowSums(Ab_counts > 0) >= dim(tumor_metadata)[1]*0.1
Ab_counts=Ab_counts[tobetaken,]
dim(Ab_counts) #484 358


#function for performing surv analysis
perform_survival_analysis <- function(mm_check, column_name) {
  test_surv_df <- data.frame(Run = mm_check$Run, Event = mm_check$RFS_Status, Time = mm_check$RFS_time_Days)
  mm_check_subset <- mm_check[, c('Run', column_name)]
  mm_check_subset <- mm_check_subset %>% arrange(desc(!!sym(column_name)))
  mm_check_subset <- mm_check_subset %>%
    mutate(Binary = ifelse(row_number() <= round(dim(mm_check_subset)[1]*0.1), "high", "low"))
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
dim(prom_df)
# write.csv(prom_df, 'all_prom.csv')
prom_sig_result <- filter(prom_df, p_value < 0.05)
dim(prom_sig_result) 
# write.csv(prom_sig_result, "prom_sig_pvalue.csv")
prom_sig_result_adj <- filter(prom_df, fdr_value < 0.1)
dim(prom_sig_result_adj) #0  3
# write.csv(prom_sig_result_adj, "prom_sig_fdr_0.1.csv")
#now doing the same thing with gene
fdr_p_values = p.adjust(unlist(gene_p_values), method = 'fdr')  
gene_df <- data.frame(sample_ID = names(gene_p_values), p_value = unlist(gene_p_values), fdr_values = fdr_p_values)
dim(gene_df)
# write.csv(gene_df, 'all_gene.csv')
gene_sig_result <- filter(gene_df, p_value < 0.05)
dim(gene_sig_result) # 39   2
write.csv(gene_sig_result, 'gene_sig_pvalue.csv')
gene_sig_result_adj <- filter(gene_df, fdr_values < 0.1)
dim(gene_sig_result_adj) #0  3
#write.csv(gene_sig_result_adj, 'gene_sig_fdr_0.1.csv')
gene_nonsig_result <- filter(gene_df, p_value > 0.05)
dim(gene_nonsig_result) # 1351    3
write.csv(gene_nonsig_result, 'gene_nonsig_pvalue.csv')
gene_nonsig_result_adj <- filter(gene_df, fdr_values > 0.1)
dim(gene_nonsig_result_adj) #1390  3
write.csv(gene_nonsig_result_adj, 'gene_nonsig_fdr_0.1.csv')

#seeing the intersection
list1 <- read_csv('prom_sig_pvalue.csv')
head(list1)
dim(list1)
list1=list1 %>% separate(sample_ID,into =c("Promoter", "Gene_names2"),sep = "_")
head(list1)

list2 = read_csv('all_gene.csv')
dim(list2)
head(list2)
length(intersect(list1$Gene_names2,list2$sample_ID))  # 1748

write_out <- intersect(list1$Gene_names2,list2$sample_ID)
length(write_out)
write.csv(write_out, 'com_p_sig_g_nonsig_fdr.csv')
# write.csv(write_out, 'com_p_sig_pv_g_nonsig_pv_fdr.csv')
write.csv(write_out, 'com_p_g_sig_pv.csv')
write.csv(write_out, 'com_p_sig_g_nonsig_pv.csv')
write.csv(write_out, 'com_p_sig_pv_g_nonsig_fdr.csv')

#for preparing the dataframe
test_count <- p_counts[,]
rownames(test_count) = test_count$ID
test_count = t(test_count)
test_count = test_count[-c(1,2,3),]
dim(surv_metadata)
test_count = test_count[rownames(test_count)%in% surv_metadata$Run , ]
dim(test_count)
test_count_2 = t(Rl_counts)
test_count_2 = test_count_2[rownames(test_count_2)%in% surv_metadata$Run , ]
dim(test_count_2)
test_count_3 = t(Ab_counts)
#let's understand the expression distribution for gene - ENSG00000136044.12 - pr12026_ENSG00000136044.12 - pr12028_ENSG00000136044.12
test_df = data.frame(gene_exp = gene_matrix$ENSG00000136044.12, prom_count = test_count[,7723], Abs_count = test_count_3[,217],
                     Rl_count =  test_count_2[,364])

#making all curves in a single frame - doesn't give clear picture
library(ggplot2)
x_limits <- range(c(test_df$gene_exp, test_df$prom_count, test_df$Abs_count, test_df$Rl_count))

ggplot(test_df, aes(x = gene_exp)) +
  geom_density(aes(color = "gene_exp"), fill = "lightblue", alpha = 0.5) +
  geom_density(aes(x = prom_count, color = "prom_count"), fill = "lightgreen", alpha = 0.5) +
  geom_density(aes(x = Abs_count, color = "Abs_count"), fill = "salmon", alpha = 0.5) +
  geom_density(aes(x = Rl_count, color = "Rl_count"), fill = "purple", alpha = 0.5) +
  scale_color_manual(name = "Expressions",
                     values = c("gene_exp"= "blue", "prom_count"= "green", "Abs_count"= "red", "Rl_count"= "purple")) +
  theme_minimal() +
  xlim(x_limits) 
  
?geom_density


# Plotting with ggplot using facets
library(ggplot2)
library(tidyr)
test_df_long <- test_df %>%
  pivot_longer(cols = c(gene_exp, prom_count, Abs_count, Rl_count),
               names_to = "variable",
               values_to = "value")
ggplot(test_df_long, aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(name = "Expressions",
                     values = c("gene_exp" = "blue", "prom_count" = "green", "Abs_count" = "red", "Rl_count" = "purple")) +
  scale_fill_manual(name = "Expressions",
                    values = c("gene_exp" = "lightblue", "prom_count" = "lightgreen", "Abs_count" = "salmon", "Rl_count" = "purple")) +
  facet_wrap(~ variable, scales = "free") +  # Create separate facets for each variable
  theme_minimal()
