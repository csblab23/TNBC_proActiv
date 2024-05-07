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
    num_ids_gt_threshold <- sum(gene_counts > 5)
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
threshold <- 0.1 * n_samples
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
dim(Ab_counts) #776 358
dim(vMat_2) #388 358

#function for performing matrix processing
perform_analysis_for_matrix <- function(metadata, matrix) {
  mm <- merge(metadata, matrix, by = 0)
  rownames(mm) <- mm$Row.names
  mm <- mm[, -1]
}

prom_p_values <- list()
gene_p_values <- list()

prom_matrix <- perform_analysis_for_matrix(surv_metadata, t(Ab_counts))
dim(prom_matrix) #358 784

gene_matrix <- perform_analysis_for_matrix(surv_metadata, t(vMat_2))
dim(gene_matrix) #358 396

gene_names = c("ENSG00000112245", "ENSG00000196704")

gene_idx = list()
sample_col_list = c()
for (gene in gene_names[1:length(gene_names)]) {
  gene_idx[[gene]] = grep(gene, filtered_df$geneId)
  for (value in gene_idx[[gene]]){
    sample_col_list = c(sample_col_list, as.character(filtered_df[value, 3]))
  }
}


for (sample_col in sample_col_list) {
  mm_check = prom_matrix[,]
  cutoff_values = c()
  mod_maxstat <- maxstat.test(Surv(RFS_time_Days, RFS_Status) ~ mm_check[, sample_col], 
                              data = mm_check, smethod = "LogRank")
  cutoff_values <- mod_maxstat$estimate
  mm_check_binary <- mm_check[, sample_col]
  mm_check_binary <- ifelse(mm_check[, sample_col] > cutoff_values[1], "high", "low")
  mm_meta <- mm_check[, 1:9]
  mm_finaldf <- mm_meta %>% mutate(sample_of_interest = mm_check_binary)
  mm_finaldf[,"sample_of_interest"] <- as.factor(mm_finaldf[,"sample_of_interest"])
  sfit <- survfit(Surv(mm_finaldf$RFS_time_Days, mm_finaldf$RFS_Status)~ mm_finaldf[,"sample_of_interest"], data=mm_finaldf)
  plot_title <- paste0("Survival Plot for ", sample_col)
  p = ggsurvplot(sfit, pval=TRUE, title= plot_title, legend.labs=c("HIGH", "LOW"), legend.title="")
  plot_object <- p$plot
  filename <- paste0("plot_", gsub(" ", "_", sample_col), ".pdf")
  ggsave(filename, plot = plot_object, device = "pdf")
}

#for plottimg gene_expression
colnames(gene_matrix)[9:dim(gene_matrix)[2]] = gsub("\\..*", "", colnames(gene_matrix)[9:dim(gene_matrix)[2]])
for (sample_col in gene_names) {
  mm_check = gene_matrix[,]
  cutoff_values = c()
  mod_maxstat <- maxstat.test(Surv(RFS_time_Days, RFS_Status) ~ mm_check[, sample_col], 
                              data = mm_check, smethod = "LogRank")
  cutoff_values <- mod_maxstat$estimate
  mm_check_binary <- mm_check[, sample_col]
  mm_check_binary <- ifelse(mm_check[, sample_col] > cutoff_values[1], "high", "low")
  mm_meta <- mm_check[, 1:9]
  mm_finaldf <- mm_meta %>% mutate(sample_of_interest = mm_check_binary)
  mm_finaldf[,"sample_of_interest"] <- as.factor(mm_finaldf[,"sample_of_interest"])
  sfit <- survfit(Surv(mm_finaldf$RFS_time_Days, mm_finaldf$RFS_Status)~ mm_finaldf[,"sample_of_interest"], data=mm_finaldf)
  plot_title <- paste0("Survival Plot for ", sample_col)
  p = ggsurvplot(sfit, pval=TRUE, title= plot_title, legend.labs=c("HIGH", "LOW"), legend.title="", palette = c("magenta", "royalblue"))
  plot_object <- p$plot
  filename <- paste0("plot_", gsub(" ", "_", sample_col), ".pdf")
  ggsave(filename, plot = plot_object, device = "pdf")
}
