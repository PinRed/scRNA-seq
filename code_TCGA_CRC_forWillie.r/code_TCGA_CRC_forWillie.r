
##################################################################
#############################s############### COAD ############################
##################################################################

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# 2. 安裝 Bioconductor 套件：limma 和 edgeR
# TCGAbiolinks 和 DESeq2 是我們前面討論過的，也建議安裝
#BiocManager::install(c("limma", "edgeR", "DESeq2", "TCGAbiolinks"), ask = FALSE)

# 3. 安裝 CRAN 套件：Seurat, matrixStats, dplyr
# Seurat 雖然主要是單細胞套件，但在某些批量分析的可視化中也可能用到。
#install.packages(c("Seurat", "matrixStats", "dplyr"), dependencies = TRUE, ask = FALSE)



<<<<<<< HEAD
#Raw.Data.Path<-"E:/Dropbox/Work_temp/IO_nanoStrings_CRC/TCGA"
#setwd(Raw.Data.Path)

COAD.Clinical<-read.table(file="survival_COAD_survival.txt",sep="\t",head=TRUE,quote = "\"")
head(COAD.Clinical)

COAD.Clinical2<-read.table(file="TCGA.COAD.sampleMap_COAD_clinicalMatrix",sep="\t",head=TRUE,quote = "\"")
head(COAD.Clinical2[,1:10])

TCGA.temp2<-read.table(file="TCGA.COAD.sampleMap_HiSeqV2.data",sep="\t",head=TRUE,quote = "\"")
head(TCGA.temp2[,1:10])
TCGA.temp<-read.table(file="TCGA.COAD.sampleMap_GAV2.data",sep="\t",head=TRUE,quote = "\"")
=======
# Raw.Data.Path<-"E:/Dropbox/Work_temp/IO_nanoStrings_CRC/TCGA"
# setwd(Raw.Data.Path)

##### View survival_COAD_survival.txt  #####
COAD.Clinical<-read.table(file="data/survival_COAD_survival.txt",sep="\t",head=TRUE,quote = "\"")
dim(COAD.Clinical)
colnames(COAD.Clinical)
head(COAD.Clinical)
colSums(is.na(COAD.Clinical[colSums(is.na(COAD.Clinical)) > 0]))



##### View TCGA.COAD.sampleMap_COAD_clinicalMatrix  #####
COAD.Clinical2<-read.table(file="data/TCGA.COAD.sampleMap_COAD_clinicalMatrix",sep="\t",head=TRUE,quote = "\"")
dim(COAD.Clinical2)
colnames(COAD.Clinical2)
head(COAD.Clinical2[,1:10])
colSums(is.na(COAD.Clinical2[colSums(is.na(COAD.Clinical2)) > 0]))


# Column categories

# 1. 基本病人資料
basic_info <- c(
  "sampleID", "patient_id", "bcr_patient_barcode", "bcr_sample_barcode", "bcr_followup_barcode",
  "gender", "age_at_initial_pathologic_diagnosis", "days_to_birth", "days_to_death",
  "days_to_last_followup", "vital_status", "height", "weight", "initial_weight",
  "year_of_initial_pathologic_diagnosis", "system_version", "form_completion_date"
)

# 2. 分子與基因相關資訊
molecular_gene <- c(
  "MSI_updated_Oct62011", "microsatellite_instability", "CIMP",
  "braf_gene_analysis_performed", "braf_gene_analysis_result",
  "kras_gene_analysis_performed", "kras_mutation_found", "kras_mutation_codon",
  "hypermutation", "loss_expression_of_mismatch_repair_proteins_by_ihc",
  "loss_expression_of_mismatch_repair_proteins_by_ihc_result",
  "non_silent_mutation", "silent_mutation", "total_mutation",
  "non_silent_rate_per_Mb", "silent_rate_per_Mb",
  "number_of_abnormal_loci", "number_of_loci_tested"
)

# 3. 臨床診斷與病理
clinical_pathology <- c(
  "primary_site", "tumor_tissue_site", "anatomic_neoplasm_subdivision",
  "histological_type", "pathologic_T", "pathologic_N", "pathologic_M",
  "pathologic_stage", "residual_tumor",
  "residual_disease_post_new_tumor_event_margin_status",
  "circumferential_resection_margin", "perineural_invasion_present",
  "lymphatic_invasion", "venous_invasion",
  "lymph_node_examined_count", "number_of_lymphnodes_positive_by_he",
  "number_of_lymphnodes_positive_by_ihc", "non_nodal_tumor_deposits", "pathology_report_file_name", "primary_lymph_node_presentation_assessment"
)

# 4. 治療與追蹤
treatment_followup <- c(
  "radiation_therapy", "additional_radiation_therapy", "additional_pharmaceutical_therapy",
  "postoperative_rx_tx", "preoperative_pretreatment_cea_level",
  "primary_therapy_outcome_success", "followup_treatment_success",
  "person_neoplasm_cancer_status",
  "new_tumor_event_after_initial_treatment", "new_neoplasm_event_type",
  "new_tumor_event_additional_surgery_procedure",
  "days_to_new_tumor_event_after_initial_treatment",
  "days_to_new_tumor_event_additional_surgery_procedure",
  "history_of_neoadjuvant_treatment", "lost_follow_up",
  "followup_case_report_form_submission_reason", "site_of_additional_surgery_new_tumor_event_mets"
)

# 5. 其他病史與檢查
history_other <- c(
  "history_of_colon_polyps", "colon_polyps_present", "synchronous_colon_cancer_present",
  "number_of_first_degree_relatives_with_cancer_diagnosis",
  "icd_10", "icd_o_3_histology", "icd_o_3_site",
  "informed_consent_verified",
  "AWG_MLH1_silencing", "AWG_cancer_type_Oct62011", "CDE_ID_3226963", "disease_code"
)

# 6. TCGA 整合與多組學欄位
tcga_multiomics <- c(
  grep("^X_", colnames(COAD.Clinical2), value = TRUE),  # 所有 X_ 開頭的欄位
  "project_code", "X_cohort", "X_primary_disease", "X_primary_site"
)

# 7. 樣本/檢體資訊
sample_info <- c(
  "days_to_collection", "days_to_initial_pathologic_diagnosis", "intermediate_dimension", "longest_dimension", "shortest_dimension",
  "is_ffpe", "oct_embedded", "sample_type", "sample_type_id",
  "tissue_prospective_collection_indicator", "tissue_retrospective_collection_indicator",
  "tissue_source_site", "vial_number"
)


# === Split to 6 data.frame ===
#COAD_basic      <- COAD.Clinical2[, basic_info]
#COAD_molecular  <- COAD.Clinical2[, molecular_gene]
#COAD_pathology  <- COAD.Clinical2[, clinical_pathology]
#COAD_treatment  <- COAD.Clinical2[, treatment_followup]
#COAD_history    <- COAD.Clinical2[, history_other]
#COAD_multiomics <- COAD.Clinical2[, tcga_multiomics]




#### TCGA.COAD.sampleMap_HiSeqV2 AND TCGA.COAD.sampleMap_GAV2 ####
<<<<<<< HEAD
TCGA.temp2<-read.table(file="data/TCGA.COAD.sampleMap_HiSeqV2.data",sep="\t",head=TRUE,quote = "\"")
head(TCGA.temp2[,1:10])
TCGA.temp<-read.table(file="data/TCGA.COAD.sampleMap_GAV2.data",sep="\t",head=TRUE,quote = "\"")
>>>>>>> 83a24d3 (Version 3)
head(TCGA.temp[,1:10])
dim(TCGA.temp)
dim(TCGA.temp2)

<<<<<<< HEAD
# TCGA.temp2[,-1]<-head(TCGA.temp2[,1:10])
# match(colnames(TCGA.temp),"TCGA.AG.3611.01")

# matrix.temp<-2^TCGA.temp[,-1]-1
# TCGA.temp2<-log2(matrix.temp*10^6/t(matrix(colSums(matrix.temp),512,60483))+1)
# head(TCGA.temp2[,1:10])

COAD.T.index<-match(unique(paste0(COAD.Clinical$X_PATIENT,"-01")),COAD.Clinical$sample,nomatch=0)
sum(COAD.T.index!=0)
COAD.N.index<-match(unique(paste0(COAD.Clinical$X_PATIENT,"-11")),COAD.Clinical$sample,nomatch=0)
sum(COAD.N.index!=0)

# COAD.T.index<-match(paste0(colnames(TCGA.temp2),"-01"),COAD.Clinical$sample,nomatch=0)
# sum(COAD.T.index!=0)

COAD.matchT.index<-match(unique(paste0(COAD.Clinical$X_PATIENT[COAD.N.index],"-01")),COAD.Clinical$sample,nomatch=0)
cbind(COAD.Clinical$sample[COAD.matchT.index],COAD.Clinical$sample[COAD.N.index])

COAD.Clinical.T<-COAD.Clinical[COAD.T.index,]

common.names<-intersect(substr(colnames(TCGA.temp),1,15),make.names(COAD.Clinical.T$sample))
common.names
length(common.names)

common.names2<-intersect(substr(colnames(TCGA.temp2),1,15),make.names(COAD.Clinical.T$sample))
common.names2
length(common.names2)

common.namesN<-intersect(substr(colnames(TCGA.temp2),1,15),make.names(COAD.Clinical$sample[COAD.N.index]))
common.namesN
length(common.namesN)


=======
=======

# Load expression data: Assumed to be log2(RSEM + 1)
coad_hiseq_logrsem<-read.table(file="data/TCGA.COAD.sampleMap_HiSeqV2.data",sep="\t",head=TRUE,quote = "\"")
head(coad_hiseq_logrsem[,1:10])
coad_gav2_logrsem<-read.table(file="data/TCGA.COAD.sampleMap_GAV2.data",sep="\t",head=TRUE,quote = "\"")
head(coad_gav2_logrsem[,1:10])
dim(coad_gav2_logrsem)
dim(coad_hiseq_logrsem)
>>>>>>> 9530aed (Add week 4)

# Note: coad_hiseq_logrsem[,1] and coad_gav2_logrsem[,1] are assumed to be Gene Names/IDs

# Match tumor cells by patient (Sample type "-01")
COAD.T.index <- match(unique(paste0(COAD.Clinical$X_PATIENT, "-01")), COAD.Clinical$sample, nomatch = 0)
sum(COAD.T.index != 0) # Count matched tumor samples

# Match normal cells by patient (Sample type "-11")
COAD.N.index <- match(unique(paste0(COAD.Clinical$X_PATIENT, "-11")), COAD.Clinical$sample, nomatch = 0)
sum(COAD.N.index != 0) # Count matched normal samples

# Find patients who have both tumor and normal cells
COAD.matchT.index <- match(unique(paste0(COAD.Clinical$X_PATIENT[COAD.N.index], "-01")), COAD.Clinical$sample, nomatch = 0)
# Combine matched tumor-normal sample pairs
cbind(COAD.Clinical$sample[COAD.matchT.index], COAD.Clinical$sample[COAD.N.index])

# Extract clinical data for tumor patients
COAD.Clinical.T <- COAD.Clinical[COAD.T.index,]

# Identify matched sample names for expression data
# 1. GAV2 Tumor samples
ga_tumor_names <- intersect(substr(colnames(coad_gav2_logrsem), 1, 15), make.names(COAD.Clinical.T$sample))
length(ga_tumor_names) 

# 2. HiSeqV2 Tumor samples
hiseq_tumor_names <- intersect(substr(colnames(coad_hiseq_logrsem), 1, 15), make.names(COAD.Clinical.T$sample))
length(hiseq_tumor_names) 

# 3. HiSeqV2 Normal samples
hiseq_normal_names <- intersect(substr(colnames(coad_hiseq_logrsem), 1, 15), make.names(COAD.Clinical$sample[COAD.N.index]))
length(hiseq_normal_names) 


#### Libraries ####
library(limma)
library(Seurat)
library(sparsepca)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

#### Data Preparation ####
# Extract expression subsets
expr_tumor_hiseq <- coad_hiseq_logrsem[, hiseq_tumor_names]
expr_normal_hiseq <- coad_hiseq_logrsem[, hiseq_normal_names]
expr_tumor_ga <- coad_gav2_logrsem[, ga_tumor_names]

# Combine HiSeqV2 Tumor and Normal data (for DEA)
expr_hiseq_tn <- cbind(expr_tumor_hiseq, expr_normal_hiseq)# Genes x Samples
group_hiseq_tn <- factor(c(rep("Tumor", ncol(expr_tumor_hiseq)), rep("Normal", ncol(expr_normal_hiseq))))

# Combine all samples data (for PCA visualization)
expr_all_samples <- cbind(expr_tumor_hiseq, expr_normal_hiseq, expr_tumor_ga)
group_all_samples <- factor(c(rep("HiSeq_Tumor", ncol(expr_tumor_hiseq)),
                              rep("HiSeq_Normal", ncol(expr_normal_hiseq)),
                              rep("GAV2_Tumor", ncol(expr_tumor_ga))))

# Filter out genes with zero variance (prevent PCA error)
expr_hiseq_filtered <- expr_hiseq_tn[rowVars(as.matrix(expr_hiseq_tn)) > 0, ]


#### PCA - All Filtered Genes ####
# Standard PCA on HiSeqV2 TN data (all filtered genes)
pca_all_genes <- prcomp(t(expr_hiseq_filtered), scale. = TRUE)

# Project all samples onto the PCA space
expr_all_filtered <- expr_all_samples[rownames(expr_hiseq_filtered), ]
expr_all_scaled <- scale(t(expr_all_filtered),
                         center = pca_all_genes$center,
                         scale = pca_all_genes$scale)
pca_all_scores <- expr_all_scaled %*% pca_all_genes$rotation

# Plot PCA
plot(pca_all_scores[,1], pca_all_scores[,2],
     col = as.numeric(group_all_samples), pch = 16,
     xlab = "PC1", ylab = "PC2", main = "PCA: All Filtered Genes")
legend("topright", legend = levels(group_all_samples),
       col = 1:length(levels(group_all_samples)), pch = 16)


#### DEA (Differential Expression Analysis) - All Filtered Genes ####

# Ensure expr_hiseq_tn is a matrix
expr_m <- as.matrix(expr_hiseq_tn)

# Create design matrix: Tumor vs Normal
design <- model.matrix(~ group_hiseq_tn)
colnames(design) <- c("Intercept", "TumorVsNormal")

# Fit linear model (limma)
fit <- limma::lmFit(expr_m, design)

# Empirical Bayes (suitable for log2 expression values)
fit <- limma::eBayes(fit, trend=TRUE)

# Get DEA results
dea_all_genes_results <- limma::topTable(fit, coef="TumorVsNormal", number=Inf, adjust.method="BH")

# Define significantly differentially expressed genes
alpha <- 0.05
logFC_cut <- 1
dea_all_genes_results$significant <- ifelse(dea_all_genes_results$adj.P.Val < alpha & abs(dea_all_genes_results$logFC) > logFC_cut, "yes", "no")
# Map gene names (from the first column of coad_hiseq_logrsem)
dea_all_genes_results$gene_name <- coad_hiseq_logrsem[rownames(dea_all_genes_results), 1]
# Subset significant genes
sig_dea_genes <- dea_all_genes_results[dea_all_genes_results$significant == "yes", ]

# Volcano plot
ggplot(dea_all_genes_results, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Tumor vs Normal (All Filtered Genes)",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  geom_text_repel(data = subset(dea_all_genes_results, significant == "yes"),
                  aes(label = gene_name), size = 3)


#### PCA - Significant DEA Genes ####
# Perform PCA using only the significant DEA genes
pca_sig_dea_genes <- prcomp(t(expr_hiseq_filtered[rownames(sig_dea_genes),]), scale. = TRUE)

# Project all samples onto the new PCA space
expr_all_dea_filtered <- expr_all_samples[rownames(sig_dea_genes), ]
expr_all_dea_scaled <- scale(t(expr_all_dea_filtered),
                             center = pca_sig_dea_genes$center,
                             scale = pca_sig_dea_genes$scale)
pca_sig_dea_scores <- expr_all_dea_scaled %*% pca_sig_dea_genes$rotation

# Plot PCA
plot(pca_sig_dea_scores[,1], pca_sig_dea_scores[,2],
     col = as.numeric(group_all_samples), pch = 16,
     xlab = "PC1", ylab = "PC2", main = "PCA: Significant DEA Genes")
legend("topright", legend = levels(group_all_samples),
       col = 1:length(levels(group_all_samples)), pch = 16)


#### Cosine Similarity: All Genes PCA vs DEA Genes PCA Loadings ####
cossim <- function(A,B) { (sum(A*B))/sqrt((sum(A^2))*(sum(B^2))) }

# Get common genes between the two PCA rotations
common_genes <- intersect(rownames(pca_all_genes$rotation), rownames(pca_sig_dea_genes$rotation))

# Extract PC loadings
rot_pca_all<- pca_all_genes$rotation[common_genes, 1:2]
rot_pca_sig_dea<- pca_sig_dea_genes$rotation[common_genes, 1:2]

# Calculate cosine similarity for PC1 and PC2 loadings
pc1_all_cos <- cossim(rot_pca_all[,1], rot_pca_sig_dea[,1])
pc2_all_cos <- cossim(rot_pca_all[,2], rot_pca_sig_dea[,2])
# cat(paste0("All vs DEA PC1 Cosine Similarity: ", pc1_all_cos, "\n"))
# cat(paste0("All vs DEA PC2 Cosine Similarity: ", pc2_all_cos, "\n"))


#### Seurat HVG Selection (500 genes) ####
# Create Seurat object (using log2(RSEM+1) for feature selection)
seurat_hiseq_tn <- CreateSeuratObject(counts = expr_hiseq_tn)
# Find 500 Highly Variable Genes (HVG)
seurat_hiseq_tn <- FindVariableFeatures(seurat_hiseq_tn, selection.method = "vst", nfeatures = 500)
hvg_genes_500 <- VariableFeatures(seurat_hiseq_tn)
expr_hvg_500 <- expr_hiseq_tn[hvg_genes_500, ]

# Perform PCA on HVGs
pca_hvg_500 <- prcomp(t(expr_hvg_500), scale. = TRUE)

# Project all samples onto the HVG PCA space
expr_hvg_all <- expr_all_samples[hvg_genes_500, ]
expr_hvg_scaled <- scale(t(expr_hvg_all),
                         center = pca_hvg_500$center,
                         scale = pca_hvg_500$scale)
pca_hvg_500_scores <- expr_hvg_scaled %*% pca_hvg_500$rotation

# Plot HVG PCA
plot(pca_hvg_500_scores[,1], pca_hvg_500_scores[,2],
     col = as.numeric(group_all_samples), pch = 16,
     xlab = "PC1", ylab = "PC2", main = "PCA: 500 HVGs")
legend("topright", legend = levels(group_all_samples),
       col = 1:length(levels(group_all_samples)), pch = 16)

#### DEA for HVGs (500 genes) ####

# Ensure expr_hvg_500 is a matrix
expr_hvg_m <- as.matrix(expr_hvg_500)

# Design matrix
design <- model.matrix(~ group_hiseq_tn)
colnames(design) <- c("Intercept", "TumorVsNormal")

# Fit linear model + Empirical Bayes
fit <- limma::lmFit(expr_hvg_m, design)
fit <- limma::eBayes(fit, trend=TRUE)

# Get DEA results
dea_hvg_500_results <- limma::topTable(fit, coef="TumorVsNormal", number=Inf, adjust.method="BH")

# Define significant genes
alpha <- 0.05
logFC_cut <- 1
dea_hvg_500_results$significant <- ifelse(dea_hvg_500_results$adj.P.Val < alpha & abs(dea_hvg_500_results$logFC) > logFC_cut, "yes", "no")

# Map gene names
dea_hvg_500_results$gene_name <- coad_hiseq_logrsem[rownames(dea_hvg_500_results), 1]

# Volcano plot
ggplot(dea_hvg_500_results, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Tumor vs Normal (HVG 500, limma DEA)",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  geom_text_repel(data = subset(dea_hvg_500_results, significant == "yes"),
                  aes(label = gene_name), size = 3)


#### PCA - Significant DEA Genes from All Genes (Re-plotted for comparison) ####
# Using the subset of significant DEA genes from All Genes PCA
pca_hvg_sig_dea <- prcomp(t(expr_hiseq_filtered[rownames(sig_dea_genes),]), scale. = TRUE)

# Project all samples onto the PCA space
expr_hvg_dea_filtered <- expr_all_samples[rownames(sig_dea_genes), ]
expr_hvg_dea_scaled <- scale(t(expr_hvg_dea_filtered),
                             center = pca_hvg_sig_dea$center,
                             scale = pca_hvg_sig_dea$scale)
pca_hvg_sig_dea_scores <- expr_hvg_dea_scaled %*% pca_hvg_sig_dea$rotation

# Plot PCA
plot(pca_hvg_sig_dea_scores[,1], pca_hvg_sig_dea_scores[,2],
     col = as.numeric(group_all_samples), pch = 16,
     xlab = "PC1", ylab = "PC2", main = "PCA: All Genes DEA (Re-plotted)")
legend("topright", legend = levels(group_all_samples),
       col = 1:length(levels(group_all_samples)), pch = 16)


#### Cosine Similarity: HVG PCA vs All Genes DEA PCA Loadings ####
# Compare HVG PCA (pca_hvg_500) and All Genes DEA PCA (pca_sig_dea_genes)
common_genes_hvg_dea <- intersect(rownames(pca_hvg_500$rotation), rownames(pca_sig_dea_genes$rotation))

# Extract PC loadings
rot_hvg_500<- pca_hvg_500$rotation[common_genes_hvg_dea, 1:2]
rot_pca_sig_dea_2<- pca_sig_dea_genes$rotation[common_genes_hvg_dea, 1:2]

# Calculate cosine similarity
pc1_hvg_dea_cos <- cossim(rot_hvg_500[,1], rot_pca_sig_dea_2[,1])
pc2_hvg_dea_cos <- cossim(rot_hvg_500[,2], rot_pca_sig_dea_2[,2])
# cat(paste0("HVG vs DEA PC1 Cosine Similarity: ", pc1_hvg_dea_cos, "\n"))
# cat(paste0("HVG vs DEA PC2 Cosine Similarity: ", pc2_hvg_dea_cos, "\n"))


#### Sparse PCA - All Filtered Genes ####
# Perform Sparse PCA on all filtered genes
spca_all_genes_result <- sparsepca::spca(t(expr_hiseq_filtered), k = 2, alpha = 0.01, max_iter = 200, scale = TRUE)

# Project all samples onto the SPCA space
expr_all_spca_filtered <- expr_all_samples[rownames(expr_hiseq_filtered), ] 
pca_spca_all_scores <- t(expr_all_spca_filtered) %*% spca_all_genes_result$loadings[, 1:2]
rownames(spca_all_genes_result$loadings) <- rownames(expr_hiseq_filtered)

# Plot Sparse PCA
pc1_loadings_count <- sum(spca_all_genes_result$loadings[,1] != 0)
pc2_loadings_count <- sum(spca_all_genes_result$loadings[,2] != 0)
plot(pca_spca_all_scores[,1], pca_spca_all_scores[,2],
     col = as.numeric(group_all_samples), pch = 16,
     xlab = "Sparse PC1", ylab = "Sparse PC2",
     main = paste0("Sparse PCA: [PC1 Loadings:", pc1_loadings_count, ", PC2 Loadings:", pc2_loadings_count, "]"))
legend("topright", legend = levels(group_all_samples),
       col = 1:length(levels(group_all_samples)), pch = 16)

# Count non-zero loadings
cat("PC1 non-zero loadings:", pc1_loadings_count, "\n")
cat("PC2 non-zero loadings:", pc2_loadings_count, "\n")
# Extract genes selected by Sparse PCA (non-zero loadings)
expr_spca_genes <- expr_hiseq_filtered[rowSums(spca_all_genes_result$loadings[,1:2] != 0) > 0, ]


#### DEA for Sparse PCA Selected Genes ####

# Ensure expr_spca_genes is a matrix
expr_spca_m <- as.matrix(expr_spca_genes)

# Design matrix
design <- model.matrix(~ group_hiseq_tn)
colnames(design) <- c("Intercept", "TumorVsNormal")

# Fit linear model + Empirical Bayes
fit <- limma::lmFit(expr_spca_m, design)
fit <- limma::eBayes(fit, trend=TRUE)

# Get DEA results
dea_spca_genes_results <- limma::topTable(fit, coef="TumorVsNormal", number=Inf, adjust.method="BH")

# Define significant genes
alpha <- 0.05
logFC_cut <- 1
dea_spca_genes_results$significant <- ifelse(dea_spca_genes_results$adj.P.Val < alpha & abs(dea_spca_genes_results$logFC) > logFC_cut, "yes", "no")

# Map gene names
dea_spca_genes_results$gene_name <- coad_hiseq_logrsem[rownames(dea_spca_genes_results), 1]

# Volcano plot
ggplot(dea_spca_genes_results, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Tumor vs Normal (Sparse PCA Genes, limma DEA)",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  geom_text_repel(data = subset(dea_spca_genes_results, significant == "yes"),
                  aes(label = gene_name), size = 3)


#### Sparse PCA - Significant DEA Genes from SPCA ####
# Extract row indices of significant DEA genes from SPCA subset
sig_dea_spca_indices <- as.character(sort(as.numeric(rownames(dea_spca_genes_results[dea_spca_genes_results$significant == "yes", ]))))
# Perform Sparse PCA on these significant genes
spca_sig_dea_result <- sparsepca::spca(t(expr_hiseq_filtered[sig_dea_spca_indices,]), k = 2, alpha = 0.01, max_iter = 200, scale = TRUE)
rownames(spca_sig_dea_result$loadings) <- sig_dea_spca_indices

# Project all samples onto the new SPCA space
pca_scores_spca <- t(expr_all_filtered[sig_dea_spca_indices,]) %*% spca_sig_dea_result$loadings[, 1:2]


# Plot SPCA
plot(pca_scores_spca[,1], pca_scores_spca[,2],
     col = as.numeric(group_all_samples), pch = 16,
     xlab = "Sparse PC1", ylab = "Sparse PC2", main = "SPCA: SPCA-DEA Significant Genes")
legend("topright", legend = levels(group_all_samples),
       col = 1:length(levels(group_all_samples)), pch = 16)

#### Cosine Similarity: All Genes SPCA vs SPCA-DEA Loadings ####
# Compare All Genes SPCA (spca_all_genes_result) and SPCA-DEA SPCA (spca_sig_dea_result)
common_genes_spca <- intersect(rownames(expr_spca_genes), sig_dea_spca_indices)

# Extract PC loadings
rot_spca_all<- spca_all_genes_result$loadings[common_genes_spca, 1:2]
rot_spca_dea<- spca_sig_dea_result$loadings[common_genes_spca, 1:2]

# Calculate cosine similarity
pc1_spca_cos <- cossim(rot_spca_all[,1], rot_spca_dea[,1])
pc2_spca_cos <- cossim(rot_spca_all[,2], rot_spca_dea[,2])
# cat(paste0("SPCA vs SPCA-DEA PC1 Cosine Similarity: ", pc1_spca_cos, "\n"))
# cat(paste0("SPCA vs SPCA-DEA PC2 Cosine Similarity: ", pc2_spca_cos, "\n"))


#### Intersection Genes Analysis ####

# Intersection of HVG (500) genes and All Genes DEA significant genes (using gene names)
intersect_hvg_dea <- intersect(coad_hiseq_logrsem[hvg_genes_500,1], sig_dea_genes$gene_name)
cat("Number of HVG (500) genes:", length(hvg_genes_500), "\n")
cat("Number of genes in HVG & DEA intersection:", length(intersect_hvg_dea), "\n")

# Intersection of Sparse PCA genes and All Genes DEA significant genes (using gene names)
spca_selected_gene_names <- coad_hiseq_logrsem[rownames(expr_spca_genes),1]
intersect_spca_dea <- intersect(spca_selected_gene_names, sig_dea_genes$gene_name)
cat("Number of Sparse PCA genes:", length(spca_selected_gene_names), "\n")
cat("Number of genes in SPCA & DEA intersection:", length(intersect_spca_dea), "\n")

# Triple intersection: HVG, SPCA, and DEA significant genes
intersect_hvg_spca_dea <- intersect(intersect_hvg_dea, intersect_spca_dea)
cat("Number of genes in Triple Intersection (HVG, SPCA, DEA):", length(intersect_hvg_spca_dea), "\n")


#### HVG (2000) + DEA + Sparse PCA Analysis ####

# 1. Find 2000 Highly Variable Genes (HVG)
seurat_hiseq_2000 <- CreateSeuratObject(counts = expr_hiseq_tn)
seurat_hiseq_2000 <- FindVariableFeatures(seurat_hiseq_2000, selection.method = "vst", nfeatures = 2000)
hvg_genes_2000 <- VariableFeatures(seurat_hiseq_2000)
expr_hvg_2000 <- expr_hiseq_tn[hvg_genes_2000, ]

# 2. Perform DEA on 2000 HVGs
expr_hvg_2000_m <- as.matrix(expr_hvg_2000)
design <- model.matrix(~ group_hiseq_tn)
colnames(design) <- c("Intercept", "TumorVsNormal")
fit <- limma::lmFit(expr_hvg_2000_m, design)
fit <- limma::eBayes(fit, trend=TRUE)
dea_hvg_2000_results <- limma::topTable(fit, coef="TumorVsNormal", number=Inf, adjust.method="BH")

# Define significant genes
alpha <- 0.05
logFC_cut <- 1
dea_hvg_2000_results$significant <- ifelse(dea_hvg_2000_results$adj.P.Val < alpha & abs(dea_hvg_2000_results$logFC) > logFC_cut, "yes", "no")
# Map gene names
dea_hvg_2000_results$gene_name <- coad_hiseq_logrsem[rownames(dea_hvg_2000_results), 1]

# Extract names of significant DEA genes among 2000 HVGs
sig_dea_hvg_2000 <- rownames(dea_hvg_2000_results[dea_hvg_2000_results$significant == 'yes',])

# 3. Perform Sparse PCA on the HVG-DEA significant genes
expr_hvg_dea_filtered_2000 <- expr_hiseq_filtered[sig_dea_hvg_2000 ,] 

spca_hvg_dea_result <- sparsepca::spca(t(expr_hvg_dea_filtered_2000), k = 2, alpha = 0.01, max_iter = 200, scale = TRUE)

# Project all samples onto the new SPCA space
expr_all_hvg_dea_filtered <- expr_all_samples[sig_dea_hvg_2000,]
pca_spca_hvg_dea_scores <- t(expr_all_hvg_dea_filtered) %*% spca_hvg_dea_result$loadings[, 1:2]
rownames(spca_hvg_dea_result$loadings) <- sig_dea_hvg_2000

# Plot Sparse PCA
pc1_loadings_count_hvgdea <- sum(spca_hvg_dea_result$loadings[,1] != 0)
pc2_loadings_count_hvgdea <- sum(spca_hvg_dea_result$loadings[,2] != 0)
plot(pca_spca_hvg_dea_scores[,1], pca_spca_hvg_dea_scores[,2],
     col = as.numeric(group_all_samples), pch = 16,
     xlab = "Sparse PC1", ylab = "Sparse PC2",
     main = paste0("HVG[2,000] + DEA [", length(sig_dea_hvg_2000), "] + Sparse PCA [", pc1_loadings_count_hvgdea, ", ", pc2_loadings_count_hvgdea, "]"))
legend("topright", legend = levels(group_all_samples),
       col = 1:length(levels(group_all_samples)), pch = 16)

# Count non-zero loadings
cat("PC1 non-zero loadings:", pc1_loadings_count_hvgdea, "\n")
cat("PC2 non-zero loadings:", pc2_loadings_count_hvgdea, "\n")
# Extract genes with non-zero loadings
spca_hvg_dea_loadings <- spca_hvg_dea_result$loadings[rowSums(spca_hvg_dea_result$loadings[,1:2] != 0) > 0, ]
spca_hvg_dea_gene_names <- coad_hiseq_logrsem[rownames(spca_hvg_dea_loadings),1]

expr_spca_hvg_dea <- expr_hiseq_filtered[rownames(spca_hvg_dea_loadings),]


#### Heatmap of Intersecting Genes ####

# Select genes that appear in the triple intersection (intersect_hvg_spca_dea)
gene_indices_to_plot <- intersect_hvg_spca_dea # Gene Names
# Get the row names (IDs) corresponding to the Gene Names
gene_indices_to_plot_rowname <- rownames(coad_hiseq_logrsem)[match(gene_indices_to_plot, coad_hiseq_logrsem[,1])]
# Subset the expression data
expr_intersect_genes <- expr_hiseq_tn[gene_indices_to_plot_rowname, ]

# Calculate the Pearson correlation matrix between genes
gene_correlation_matrix <- cor(t(expr_intersect_genes)) 

# Rename rows/columns with Gene Names
rownames(gene_correlation_matrix) <- gene_indices_to_plot
colnames(gene_correlation_matrix) <- gene_indices_to_plot

dim(gene_correlation_matrix)

# Plot heatmap
pheatmap(gene_correlation_matrix,
         main = paste0("Gene Correlation Heatmap (Intersecting ", length(intersect_hvg_spca_dea), " Genes)"),
         clustering_distance_rows = "euclidean", # Gene clustering distance
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6, # Adjust font size for large number of genes
         fontsize_col = 6
)



<<<<<<< HEAD
#### Analysis ####
>>>>>>> 83a24d3 (Version 3)
=======

### LASSO on Expression Data (Tumor vs Normal)

library(glmnet)

# expr_filtered: 基因 x 樣本 (經過過濾零方差的基因)
# group_main: 樣本分組 (Tumor / Normal)

X <- as.matrix(t(expr_hvg_2000))  # glmnet 要用樣本 x 特徵
group_main <- factor(c(rep("Tumor",  ncol(expr_tumor_hiseq)),
                       rep("Normal", ncol(expr_normal_hiseq))
                      ))
y <- as.factor(group_main)

# 轉成二元 0/1 (Tumor=1, Normal=0)
y_bin <- ifelse(y == "Tumor", 1, 0)

# 建立 LASSO logistic regression
set.seed(123)
cvfit <- cv.glmnet(X, y_bin, family="binomial", alpha=1, nfolds=10)

# 找最佳 lambda
best_lambda <- cvfit$lambda.1se
cat("Best lambda:", best_lambda, "\n")

# 用最佳 lambda 重新訓練
lasso_fit <- glmnet(X, y_bin, family="binomial", alpha=1, lambda=best_lambda)

# 提取非零係數基因 (被 LASSO 選中的特徵)
coef_lasso <- coef(lasso_fit)
selected_genes <- rownames(coef_lasso)[which(coef_lasso != 0)]
selected_genes <- selected_genes[-1]  # 去掉 intercept

cat("Number of selected genes:", length(selected_genes), "\n")
head(selected_genes)
coad_hiseq_logrsem[selected_genes,1]

>>>>>>> 9530aed (Add week 4)
