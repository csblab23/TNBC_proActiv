setwd("~/univariate_surv_proActiv_new/univariate_surv_proActiv/loop_v2_apr16")

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
result2 <- readRDS('result2_v2_updated.rds')
metadata <- as.data.frame(read_csv("../meta_with_batchinfo.csv"))
count <- as.data.frame(read_csv("../Featurecounts_combined_fuscc.csv"))

# Data preprocessing
metadata[1:5,1:5]
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run

tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor')
dim(tumor_metadata) # 360  99
table(tumor_metadata$Tumor_Normal1)

# taking >=2 promoters gene
result2_promoter = rowData(result2)
filtered_df <- result2_promoter %>%
  as.data.frame() %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep=''))

filtered_df <- filtered_df[, c('promoterId','geneId', 'ID')]
dim(filtered_df) #24120     3
length(unique(filtered_df$geneId)) # 9315

table(table(filtered_df$geneId))

promoter_counts <- as.data.frame(result2@assays@data@listData$promoterCounts)
length(intersect(filtered_df$promoterId, rownames(promoter_counts)))  ## 24120
keep = intersect(filtered_df$promoterId, rownames(promoter_counts))
promoter_counts <- promoter_counts[keep,]
dim(promoter_counts) # 24120   360
#taking tumor sample only
promoter_counts[1:5,1:5]
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

# for taking genes that are present in atleast 10% of samples
n_samples <- length(gene_lists)
ensg_freq <- table(unlist(gene_lists))
threshold <- round(0.1 * n_samples)
ensg_10percent <- names(ensg_freq[ensg_freq >= threshold])
length(ensg_10percent) #  2378
dim(filtered_df) # 24120     3
length(unique(filtered_df$geneId)) # 9315
filtered_df <- filtered_df[filtered_df$geneId %in% ensg_10percent,]
dim(filtered_df) # 7162    3
length(unique(filtered_df$geneId)) # 2378
table(table(filtered_df$geneId))
#    2    3    4    5    6    7    8    9   10   11   12   16   18   33
# 1163  661  291  131   53   30   20    9    4    6    5    3    1    1



#takinng the genes that have promoter with relative activity > 0.5 in at least %10 of the tumor samples
Rl_counts <- as.data.frame(result2@assays@data@listData$relativePromoterActivity)
dim(Rl_counts) # 53910   360
Rl_counts <- Rl_counts[complete.cases(Rl_counts),]
dim(Rl_counts) #  8506  360
length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  # 3903
keep = intersect(filtered_df$promoterId, rownames(Rl_counts))
Rl_counts <- Rl_counts[keep,]
dim(Rl_counts) # 3903  360

#taking tumor sample only
Rl_counts[1:5,1:5]
colnames(Rl_counts) <- gsub("_SJ.out", "", colnames(Rl_counts))
Rl_counts <- Rl_counts[,tumor_metadata$Run]
dim(Rl_counts) # 3903  360

#TAKING only the promoters with relative activity > 0.5 in at least %10 of the tumor samples
tobetaken_rel <- rowSums(Rl_counts > 0.5) >= dim(Rl_counts)[2]*0.1
table(tobetaken_rel)
# #tobetaken_rel
# FALSE  TRUE
#  2415  1488  
Rl_counts=Rl_counts[tobetaken_rel,]
dim(Rl_counts) # 1488  360
length(intersect(filtered_df$promoterId, rownames(Rl_counts)))  ## 1488
filtered_df <- filtered_df[filtered_df$promoterId %in% rownames(Rl_counts),]
dim(filtered_df) # 1488   3
length(unique(filtered_df$geneId)) # 1099
table(table(filtered_df$geneId))

# Take into consideration genes that have more than 2 promoters after rel PA filter
ensg_freq_2 <- table(filtered_df$geneId)
ensg_final <- names(ensg_freq_2[ensg_freq_2 >= 2])
length(ensg_final) # 389
filtered_df <- filtered_df[filtered_df$geneId %in% ensg_final,]
dim(filtered_df) # 778   3
length(unique(filtered_df$geneId)) # 389
table(table(filtered_df$geneId))
#   2
# 389

# pre-processing on feature count matrix
count <- read_csv("../Featurecounts_combined_fuscc.csv") %>%
  as.data.frame() %>%
  column_to_rownames("Geneid") %>%
  subset(., select = -1) %>%
  filter(rownames(.) %in% filtered_df$geneId)
count[1:5,1:5]
dim(count) # 389 448
# remove the normal samples
count=count[,colnames(count) %in% tumor_metadata$Run]
dim(count) #  388 360
count=count[,tumor_metadata$Run]
all(tumor_metadata$Run==colnames(count))


# Apply cpm
library(edgeR)
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= dim(tumor_metadata)[1]*0.1
count <- count[is.exprs, ]
dim(count) #  388 360
x <- DGEList(counts=count )

x <- calcNormFactors(x,method = "TMM")
v <- voom(x, plot=F)
vMat <- v$E
dim(vMat)

#batch adjustment
library(sva)
vMat_2 <- sva::ComBat(vMat, batch=tumor_metadata$batch)
vMat_2 <- vMat_2[, tumor_metadata$Run] #removing normal sample
dim(vMat_2) # 388 360

p_info <- filtered_df[filtered_df$geneId %in% row.names(vMat_2),]
dim(p_info) # 776   3
length(unique(p_info$geneId)) # 388

#taking Abs PA
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
Ab_counts <- Ab_counts[rownames(Ab_counts) %in% p_info$promoterId,]
Ab_counts[1:5,1:5]
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
dim(vMat_2) #  388 358
dim(Ab_counts) # 776 358
tobetaken <- rowSums(Ab_counts > 0) >= dim(tumor_metadata)[1]*0.1
Ab_counts=Ab_counts[tobetaken,]
dim(Ab_counts) # 776 358

#function for performing matrix processing
perform_analysis_for_matrix <- function(metadata, matrix) {
  mm <- merge(metadata, matrix, by = 0)
  rownames(mm) <- mm$Row.names
  mm <- mm[, -1]
}
prom_matrix <- perform_analysis_for_matrix(surv_metadata, t(Ab_counts))
dim(prom_matrix) #358 875

prom_p_values <- list()
library(maxstat)
library(survival)

analyze_gene_expression <- function(mm_check) {
  cutoff_values <- numeric(length = ncol(mm_check[, 9:dim(mm_check)[2]]))
  for (i in 9:dim(mm_check)[2]) {
    mod_maxstat <- maxstat.test(Surv(RFS_time_Days, RFS_Status) ~ mm_check[, i],
                                data = mm_check, smethod = "LogRank", pmethod = c("none"))
    cutoff_values[i] <- mod_maxstat$estimate
    cat("cut-off values", "gene", i, ":", cutoff_values[i], "\n")
  }
  mm_check_binary <- mm_check[, 9:dim(mm_check)[2]]
  for (i in 1:ncol(mm_check_binary)) {
    mm_check_binary[, i] <- ifelse(mm_check[, i + 8] > cutoff_values[i + 8], "high", "low")
  }
  mm_meta <- mm_check[, 1:8]
  mm_finaldf <- merge(mm_meta, mm_check_binary, by = 0)
  mm_finaldf <- tibble::column_to_rownames(mm_finaldf, "Row.names")
  for (i in 9:dim(mm_finaldf)[2]) {
    if (is.character(mm_finaldf[, i])) {
      mm_finaldf[, i] <- factor(mm_finaldf[, i])
    }
  }
  gene_names <- colnames(mm_finaldf)[9:dim(mm_finaldf)[2]]
  gene_p_values1 <- numeric(length = length(gene_names))
  for (i in 1:length(gene_names)) {
    gene_expression1 <- mm_finaldf[, i + 8]
    fit <- coxph(Surv(mm_finaldf$RFS_time_Days, mm_finaldf$RFS_Status) ~ gene_expression1,
                 data = mm_finaldf)
    gene_p_values1[i] <- summary(fit)$coefficients['gene_expression1low', "Pr(>|z|)"]
    cat("P-value for", "gene", i, ":", gene_p_values1[i], "\n")
  }
  return(gene_p_values1)
}

#OUT FILES FOR PROMOTERS
prom_p_values <- analyze_gene_expression(prom_matrix)
prom_fdr_p_values = p.adjust(prom_p_values, method = 'fdr')
result_prom <- data.frame(Gene_names = colnames(prom_matrix)[c(9:dim(prom_matrix)[2])], Adj_pvalue = prom_fdr_p_values, PValue = prom_p_values)
dim(result_prom) #176   3
result_prom=result_prom %>% separate(Gene_names,into =c("Promoter", "Gene_names2"),sep = "_", remove = FALSE)
write.csv(result_prom, 'all_prom.csv')

sig_result <- filter(result_prom, PValue < 0.05)
dim(sig_result) # 89  4
length(unique(sig_result$Gene_names2)) #65
write.csv(sig_result, 'prom_sig_pvalue.csv')

sig_result_adj <- filter(result_prom, Adj_pvalue < 0.05)
dim(sig_result_adj) # 43 5
length(unique(sig_result_adj$Gene_names2)) #37
write.csv(sig_result_adj, 'prom_sig_fdr.csv')

#now performing the same thing for gene
perform_analysis_for_matrix <- function(metadata, matrix) {
  mm <- merge(metadata, matrix, by = 0)
  rownames(mm) <- mm$Row.names
  mm <- mm[, -1]
}

gene_matrix <- perform_analysis_for_matrix(surv_metadata, t(vMat_2))
dim(gene_matrix) #113  97
gene_p_values <- list()
gene_p_values <- analyze_gene_expression(gene_matrix)
gene_fdr_p_values = p.adjust(gene_p_values, method = 'fdr')
result_gene <- data.frame(Gene_names = colnames(gene_matrix)[c(9:dim(gene_matrix)[2])], Adj_pvalue = gene_fdr_p_values, PValue = gene_p_values)
dim(result_gene) #388   3

write.csv(result_gene, 'all_gene.csv')

sig_result_gene <- filter(result_gene, PValue < 0.1)
dim(sig_result_gene) # 35 4
length(unique(sig_result_gene$Gene_names)) #65
write.csv(sig_result_gene, 'gene_sig_pv.csv')

sig_result_adj_gene <- filter(result_gene, Adj_pvalue < 0.1)
dim(sig_result_adj_gene) # 34  3
write.csv(sig_result_adj_gene, 'gene_sig_fdr.csv')

nonsig_result_gene <- filter(result_gene, PValue > 0.1)
dim(nonsig_result_gene) # 89  4
length(unique(nonsig_result_gene$Gene_names)) #65
write.csv(nonsig_result_gene, 'gene_nonsig_pv.csv')

nonsig_result_adj_gene <- filter(result_gene, Adj_pvalue > 0.1)
dim(nonsig_result_adj_gene) # 79  3
write.csv(nonsig_result_adj_gene, 'gene_nonsig_fdr.csv')


#common genes
list1 = read_csv('prom_sig_pvalue.csv')
dim(list1)
head(list1)
length(unique(list1$Gene_names2))

list2 = read_csv('gene_nonsig_pv.csv')
dim(list2)
head(list2)

length(intersect(list1$Gene_names2, list2$Gene_names))
keep = intersect(list1$Gene_names2, list2$Gene_names)
length(keep)
keep

write.csv(keep, 'com_p_sig_g_nonsig_pv.csv')
