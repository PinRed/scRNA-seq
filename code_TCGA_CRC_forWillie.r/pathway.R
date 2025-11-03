#### COAD Pathway Analysis ####


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


#### TCGA.COAD.sampleMap_HiSeqV2 AND TCGA.COAD.sampleMap_GAV2 ####

# Load expression data: Assumed to be log2(RSEM + 1)
coad_hiseq_rsem<-read.table(file="data/TCGA.COAD.sampleMap_HiSeqV2.data",sep="\t",head=TRUE,quote = "\"")
rownames(coad_hiseq_rsem)<-coad_hiseq_rsem$sample
head(coad_hiseq_rsem[,1:10])
coad_gav2_rsem<-read.table(file="data/TCGA.COAD.sampleMap_GAV2.data",sep="\t",head=TRUE,quote = "\"")
rownames(coad_gav2_rsem)<-coad_gav2_rsem$sample
head(coad_gav2_rsem[,1:10])
dim(coad_gav2_rsem)
dim(coad_hiseq_rsem)


# Note: coad_hiseq_rsem[,1] and coad_gav2_rsem[,1] are assumed to be Gene Names/IDs

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
ga_tumor_names <- intersect(substr(colnames(coad_gav2_rsem), 1, 15), make.names(COAD.Clinical.T$sample))
length(ga_tumor_names) 

# 2. HiSeqV2 Tumor samples
hiseq_tumor_names <- intersect(substr(colnames(coad_hiseq_rsem), 1, 15), make.names(COAD.Clinical.T$sample))
length(hiseq_tumor_names) 

# 3. HiSeqV2 Normal samples
hiseq_normal_names <- intersect(substr(colnames(coad_hiseq_rsem), 1, 15), make.names(COAD.Clinical$sample[COAD.N.index]))
length(hiseq_normal_names) 



# Extract expression subsets
expr_tumor_hiseq <- coad_hiseq_rsem[, hiseq_tumor_names]
expr_normal_hiseq <- coad_hiseq_rsem[, hiseq_normal_names]
expr_tumor_ga <- coad_gav2_rsem[, ga_tumor_names]

# Loading Package
# 請先安裝 / 載入必要套件
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

# --- 檢查與合併表達矩陣 ---
# 假設 expr_normal_hiseq, expr_tumor_hiseq 已存在且為 numeric matrix / data.frame
stopifnot(is.matrix(expr_normal_hiseq) || is.data.frame(expr_normal_hiseq))
stopifnot(is.matrix(expr_tumor_hiseq) || is.data.frame(expr_tumor_hiseq))
expr_data <- cbind(expr_normal_hiseq, expr_tumor_hiseq)

# 確認 rownames (基因 ID) 存在
if(is.null(rownames(expr_data))){
  stop("expr_data lacks rownames (gene IDs). 請確認每列有基因 ID (例如 SYMBOL)。")
}

# 建立樣本群組（levels 決定設計矩陣順序）
sample_groups <- factor(c(rep("Normal", ncol(expr_normal_hiseq)), 
                          rep("Tumor", ncol(expr_tumor_hiseq))),
                        levels = c("Normal", "Tumor"))

design <- model.matrix(~0 + sample_groups)
colnames(design) <- c("Normal", "Tumor")
if(ncol(design) != 2) stop("Design matrix columns != 2，請檢查 sample_groups")

# --- 線性模型與 eBayes ---
fit <- lmFit(expr_data, design)
contrast.matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# --- 取得 topTable（取 t 作排序依據） ---
dea_results <- topTable(fit2, coef="Tumor_vs_Normal", number=Inf, sort.by="none")
# 確認 rownames 與你的 ID type
message("前 6 個 rownames（基因 ID）：")
print(head(rownames(dea_results), 6))

# 以 t 值為排名指標（也可換成 logFC）
gene_list <- dea_results$t
names(gene_list) <- rownames(dea_results)
gene_list_sorted <- sort(gene_list, decreasing = TRUE)

# --- SYMBOL -> ENTREZID 轉換（並健全性檢查） ---
# 先檢查是否為 SYMBOL（若你的 rownames 不是 SYMBOL，請改用 fromType = "ENSEMBL" 等）
unique_ids_sample <- head(names(gene_list_sorted))
message("示例基因 ID (前6)：", paste(unique_ids_sample, collapse = ", "))

s2e <- bitr(names(gene_list_sorted),
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

if(nrow(s2e) == 0) stop("bitr 未找到任何映射，請確認 rownames 是 SYMBOL 或改用正確的 fromType。")

# 將 dea_results (含 t) 與 s2e 合併，並處理重複的 ENTREZID（保留絕對 t 值最大者）
df_rank <- data.frame(SYMBOL = names(gene_list_sorted),
                      t = as.numeric(gene_list_sorted),
                      stringsAsFactors = FALSE)

df_map <- merge(df_rank, s2e, by = "SYMBOL")
# 若一個 ENTREZID 對應多個 SYMBOL，選擇絕對 t 值最大的那個來代表該 ENTREZID
df_map_unique <- df_map %>%
  group_by(ENTREZID) %>%
  arrange(desc(abs(t))) %>%
  slice(1) %>%
  ungroup()

# 建命名向量：names = ENTREZID (as character)，values = t
gene_list_entrez_sorted <- df_map_unique$t
names(gene_list_entrez_sorted) <- as.character(df_map_unique$ENTREZID)
gene_list_entrez_sorted <- sort(gene_list_entrez_sorted, decreasing = TRUE)

# 檢查長度
message("最後用於 GSEA 的基因數: ", length(gene_list_entrez_sorted))

# --- 取得 KEGG gene sets (msigdbr) 並確保 ID 型態一致 (string) ---
m_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
if(nrow(m_df) == 0) stop("msigdbr 沒抓到任何 KEGG 路徑，請確認 collection/subcollection 名稱。")
kegg_gene_sets <- m_df %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(as.character(ncbi_gene))) %>%
  { setNames(.$genes, .$gs_name) }  # 變成 named list: names = gs_name, values = character vector of ENTREZ

message("已成功獲取 KEGG 路徑數量: ", length(kegg_gene_sets))

# 篩掉基因數過少或過多的 pathway（可選）
kegg_gene_sets_filtered <- kegg_gene_sets[sapply(kegg_gene_sets, length) >= 10 & sapply(kegg_gene_sets, length) <= 500]

# --- 執行 fgsea ---
fgsea_res <- fgsea(pathways = kegg_gene_sets_filtered,
                   stats = gene_list_entrez_sorted,
                   minSize = 10,
                   maxSize = 500,
                   nperm = 1000)

# 整理結果 Sorting by FDR and NES
fgsea_res_df <- as.data.frame(fgsea_res) %>% arrange(padj,desc(abs(NES)))
message("顯示 GSEA 前 10 條結果：")
print(head(fgsea_res_df, 10))

# --- 繪圖: 若要畫單一路徑建議用 plotEnrichment 或 plotGseaTable with correct subset ---
k = 4
if(nrow(fgsea_res_df) > 0){
  top_pathway <- as.character(fgsea_res_df$pathway[k])
  message("繪製: ", top_pathway)
  # plotEnrichment expects pathway (vector of gene IDs) and stats named by gene IDs
  p <- plotEnrichment(kegg_gene_sets_filtered[[top_pathway]], gene_list_entrez_sorted) + ggtitle(top_pathway)
  print(p)
}


library(dplyr)

top5_pathways <- fgsea_res_df %>%
  arrange(padj, desc(abs(NES))) %>%
  mutate(leadingEdgeCount = sapply(leadingEdge,length)) %>%
  select(pathway, padj, NES, size, leadingEdgeCount) %>%
  head(5)
print(top5_pathways)

###### log FC #####################

# 以 logFC 為排名指標（也可換成 t 值）
gene_list_logFC <- dea_results$logFC
names(gene_list_logFC) <- rownames(dea_results)
gene_list_sorted_logFC <- sort(gene_list_logFC, decreasing = TRUE)

# --- SYMBOL -> ENTREZID 轉換（並健全性檢查） ---
# 先檢查是否為 SYMBOL（若你的 rownames 不是 SYMBOL，請改用 fromType = "ENSEMBL" 等）
unique_ids_sample_logFC <- head(names(gene_list_sorted_logFC))
message("示例基因 ID (前6)：", paste(unique_ids_sample_logFC, collapse = ", "))

s2e_logFC <- bitr(names(gene_list_sorted_logFC),
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

if(nrow(s2e) == 0) stop("bitr 未找到任何映射，請確認 rownames 是 SYMBOL 或改用正確的 fromType。")

# 將 dea_results (含 t) 與 s2e 合併，並處理重複的 ENTREZID（保留絕對 t 值最大者）
df_rank_logFC <- data.frame(SYMBOL = names(gene_list_sorted_logFC),
                      t = as.numeric(gene_list_sorted_logFC),
                      stringsAsFactors = FALSE)

df_map_logFC <- merge(df_rank_logFC, s2e_logFC, by = "SYMBOL")
# 若一個 ENTREZID 對應多個 SYMBOL，選擇絕對 t 值最大的那個來代表該 ENTREZID
df_map_unique_logFC <- df_map_logFC %>%
  group_by(ENTREZID) %>%
  arrange(desc(abs(t))) %>%
  slice(1) %>%
  ungroup()

# 建命名向量：names = ENTREZID (as character)，values = t
gene_list_entrez_sorted_logFC <- df_map_unique_logFC$t
names(gene_list_entrez_sorted_logFC) <- as.character(df_map_unique_logFC$ENTREZID)
gene_list_entrez_sorted_logFC <- sort(gene_list_entrez_sorted_logFC, decreasing = TRUE)

# 檢查長度
message("最後用於 GSEA 的基因數: ", length(gene_list_entrez_sorted_logFC))

# --- 取得 KEGG gene sets (msigdbr) 並確保 ID 型態一致 (string) ---
m_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
if(nrow(m_df) == 0) stop("msigdbr 沒抓到任何 KEGG 路徑，請確認 collection/subcollection 名稱。")
kegg_gene_sets <- m_df %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(as.character(ncbi_gene))) %>%
  { setNames(.$genes, .$gs_name) }  # 變成 named list: names = gs_name, values = character vector of ENTREZ

message("已成功獲取 KEGG 路徑數量: ", length(kegg_gene_sets))

# 篩掉基因數過少或過多的 pathway（可選）
kegg_gene_sets_filtered <- kegg_gene_sets[sapply(kegg_gene_sets, length) >= 10 & sapply(kegg_gene_sets, length) <= 500]

# --- 執行 fgsea ---
fgsea_res_logFC <- fgsea(pathways = kegg_gene_sets_filtered,
                   stats = gene_list_entrez_sorted_logFC,
                   minSize = 10,
                   maxSize = 500,
                   nperm = 1000)

# 整理結果
fgsea_res_df_logFC <- as.data.frame(fgsea_res_logFC) %>% arrange(padj,desc(abs(NES)))
message("顯示 GSEA 前 10 條結果：")
print(head(fgsea_res_df_logFC, 10))

# --- 繪圖: 若要畫單一路徑建議用 plotEnrichment 或 plotGseaTable with correct subset ---
k = 1
if(nrow(fgsea_res_df_logFC) > 0){
  top_pathway_logFC <- as.character(fgsea_res_df_logFC$pathway[k])
  message("繪製: ", top_pathway_logFC)
  # plotEnrichment expects pathway (vector of gene IDs) and stats named by gene IDs
  p <- plotEnrichment(kegg_gene_sets_filtered[[top_pathway_logFC]], gene_list_entrez_sorted_logFC) + ggtitle(paste("logFC top 1 : " ,top_pathway_logFC))
  print(p)
}


k = 1
if(nrow(fgsea_res_df) > 0){
  top_pathway <- as.character(fgsea_res_df$pathway[k])
  message("繪製: ", top_pathway)
  # plotEnrichment expects pathway (vector of gene IDs) and stats named by gene IDs
  p <- plotEnrichment(kegg_gene_sets_filtered[[top_pathway]], gene_list_entrez_sorted) + ggtitle(paste("t-value top 1 : " ,top_pathway))
  print(p)
}



library(dplyr)

top5_pathways_logFC <- fgsea_res_df_logFC %>%
  arrange(padj, desc(abs(NES))) %>%
  select(pathway, padj, NES, size, leadingEdge) %>%
  head(5)
print(top5_pathways_logFC)

##### Common pathway ####

topk_pathways_logFC <- fgsea_res_df_logFC %>%
  arrange(padj, desc(abs(NES))) %>%
  # 新增 leadingEdgeCount 欄位
  mutate(leadingEdgeCount = sapply(leadingEdge, length)) %>%
  select(pathway, padj, NES, size, leadingEdgeCount) 
print(topk_pathways_logFC)

topk_pathways <- fgsea_res_df %>%
  arrange(padj, desc(abs(NES))) %>%
  mutate(leadingEdgeCount = sapply(leadingEdge,length)) %>%
  select(pathway, padj, NES, size, leadingEdgeCount) 
print(topk_pathways)

# common leadingedge for each pathway
# 假設 common_pathways 是兩邊都有的共同 pathway
common_leadingedge_intersect <- lapply(common_pathways, function(pw) {
  le1 <- top_common_leadingedge %>% filter(pathway == pw) %>% pull(leadingEdge) %>% unlist()
  le2 <- top_common_logFC_leadingedge %>% filter(pathway == pw) %>% pull(leadingEdge) %>% unlist()
  common_genes <- intersect(le1, le2)
  
  tibble(
    pathway = pw,
    common_leadingEdge = list(common_genes),
    leadingEdgeCount = length(common_genes)
  )
}) %>% bind_rows()


# Common pathway choose padj < 0.05
p_0.05_pathways <- topk_pathways[topk_pathways$padj<0.05,]
p_0.05_pathways_logFC <- topk_pathways_logFC[topk_pathways_logFC$padj<0.05,]


common_pathways <- intersect(p_0.05_pathways$pathway, p_0.05_pathways_logFC$pathway)

# Show all columns for these overlapping pathways in p_0.05_logFC_common
p_0.05_logFC_common <- p_0.05_pathways_logFC %>%
  filter(pathway %in% common_pathways)

print(p_0.05_logFC_common)

# Similarly, for topk_pathways (T value)
p_0.05_T_common <- p_0.05_pathways %>%
  filter(pathway %in% common_pathways)

print(p_0.05_T_common)

p_0.05_common <- p_0.05_T_common
idx <- match(p_0.05_common$pathway, p_0.05_logFC_common$pathway)

# 將對應的 leadingEdgeCount 加入
p_0.05_common$leadingEdgeCount_logFC <- p_0.05_logFC_common$leadingEdgeCount[idx]

# common_pathways top 10
common_pathways <- intersect(topk_pathways_logFC[1:10,1], topk_pathways[1:10,1])

p_0.05_common_top <- p_0.05_common %>%
  filter(pathway %in% common_pathways)

# plot top 10 t value
pathway_name <- fgsea_res_df$pathway[1:k]

plotGseaTable(
  pathways = kegg_gene_sets_filtered[pathway_name],
  stats = gene_list_entrez_sorted,
  fgseaRes = fgsea_res,
  gseaParam = 0.5
)

# plot top 10 logFC
pathway_name_logFC <- fgsea_res_df_logFC$pathway[1:10]

plotGseaTable(
  pathways = kegg_gene_sets_filtered[pathway_name_logFC],
  stats = gene_list_entrez_sorted_logFC,
  fgseaRes = fgsea_res_logFC,
  gseaParam = 0.5
)




#### Random Choose ####

# --- SYMBOL -> ENTREZID 轉換（並健全性檢查） ---
# 先檢查是否為 SYMBOL（若你的 rownames 不是 SYMBOL，請改用 fromType = "ENSEMBL" 等）
unique_ids_sample <- head(names(gene_list_sorted))
message("示例基因 ID (前6)：", paste(unique_ids_sample, collapse = ", "))

s2e <- bitr(names(gene_list_sorted),
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

if(nrow(s2e) == 0) stop("bitr 未找到任何映射，請確認 rownames 是 SYMBOL 或改用正確的 fromType。")

# 將 dea_results (含 t) 與 s2e 合併，並處理重複的 ENTREZID（保留絕對 t 值最大者）
df_rank <- data.frame(SYMBOL = names(gene_list_sorted),
                      t = as.numeric(gene_list_sorted),
                      stringsAsFactors = FALSE)

df_map <- merge(df_rank, s2e, by = "SYMBOL")
# 若一個 ENTREZID 對應多個 SYMBOL，選擇絕對 t 值最大的那個來代表該 ENTREZID
df_map_unique <- df_map %>%
  group_by(ENTREZID) %>%
  arrange(desc(abs(t))) %>%
  slice(1) %>%
  ungroup()

# 建命名向量：names = ENTREZID (as character)，values = t
gene_list_entrez_sorted <- df_map_unique$t
names(gene_list_entrez_sorted) <- as.character(df_map_unique$ENTREZID)
gene_list_entrez_sorted <- sort(gene_list_entrez_sorted, decreasing = TRUE)

# 檢查長度
message("最後用於 GSEA 的基因數: ", length(gene_list_entrez_sorted))

# --- 取得 KEGG gene sets (msigdbr) 並確保 ID 型態一致 (string) ---
m_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
if(nrow(m_df) == 0) stop("msigdbr 沒抓到任何 KEGG 路徑，請確認 collection/subcollection 名稱。")
kegg_gene_sets <- m_df %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(as.character(ncbi_gene))) %>%
  { setNames(.$genes, .$gs_name) }  # 變成 named list: names = gs_name, values = character vector of ENTREZ

message("已成功獲取 KEGG 路徑數量: ", length(kegg_gene_sets))

# 篩掉基因數過少或過多的 pathway（可選）
kegg_gene_sets_filtered <- kegg_gene_sets[sapply(kegg_gene_sets, length) >= 10 & sapply(kegg_gene_sets, length) <= 500]


k = 100
r_values <- c(0.4, 0.6, 0.8, 1)
pathway_counts_all <- list()

for (r in r_values){
  message("=== Running r = ", r, " ===")
  pathway_counts <- list()
  
  for (i in 1:k){
    print(i)
    num_normal <- r*length(expr_normal_hiseq)
    num_tumor <- r*length(expr_tumor_hiseq)
    rand_normal_hiseq <- sample(expr_normal_hiseq,num_normal)
    rand_tumor_hiseq <- sample(expr_tumor_hiseq,num_tumor)
    expr_data <- cbind(rand_normal_hiseq, rand_tumor_hiseq)
    
    # 確認 rownames (基因 ID) 存在
    if(is.null(rownames(expr_data))){
      stop("expr_data lacks rownames (gene IDs). 請確認每列有基因 ID (例如 SYMBOL)。")
    }
    
    # 建立樣本群組（levels 決定設計矩陣順序）
    sample_groups <- factor(c(rep("Normal", ncol(rand_normal_hiseq)), 
                              rep("Tumor", ncol(rand_tumor_hiseq))),
                            levels = c("Normal", "Tumor"))
    
    design <- model.matrix(~0 + sample_groups)
    colnames(design) <- c("Normal", "Tumor")
    if(ncol(design) != 2) stop("Design matrix columns != 2，請檢查 sample_groups")
    
    # --- 線性模型與 eBayes ---
    fit <- lmFit(expr_data, design)
    contrast.matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # --- 取得 topTable（取 t 作排序依據） ---
    dea_results <- topTable(fit2, coef="Tumor_vs_Normal", number=Inf, sort.by="none")
    # 確認 rownames 與你的 ID type
    #message("前 6 個 rownames（基因 ID）：")
    #print(head(rownames(dea_results), 6))
    
    # 以 t 值為排名指標（也可換成 logFC）
    gene_list <- dea_results$t
    names(gene_list) <- rownames(dea_results)
    gene_list_sorted <- sort(gene_list, decreasing = TRUE)
    
    # --- 執行 fgsea ---
    fgsea_res <- fgsea(pathways = kegg_gene_sets_filtered,
                       stats = gene_list_entrez_sorted,
                       minSize = 10,
                       maxSize = 500,
                       nperm = 1000)
    
    # 整理結果 Sorting by FDR and NES
    fgsea_res_df <- as.data.frame(fgsea_res) %>% arrange(padj,desc(abs(NES)))
    
    # 找出 padj 最大的 pathway（若有 ties 會抓全部）
    max_padj_pathways <- fgsea_res_df %>%
      filter(padj == min(padj, na.rm = TRUE)) %>%
      pull(pathway) %>%
      na.omit() %>%
      unique()
    print(max_padj_pathways)
    
    # 若有 pathway 值才執行
    if (length(max_padj_pathways) > 0) {
      for (pw in max_padj_pathways) {
        if (!is.null(pw) && pw != "") {
          if (!pw %in% names(pathway_counts)) {
            pathway_counts[[pw]] <- 1
          } else {
            pathway_counts[[pw]] <- pathway_counts[[pw]] + 1
          }
        }
      }
    }
  }
  pathway_counts_all[[paste0("r_", r)]] <- pathway_counts
}



for (r in r_values) {
  r_name <- paste0("r_", r)          # 生成對應的 list 名稱
  counts <- table(unlist(pathway_counts_all[[r_name]]))
  print(paste("Table for", r_name))
  print(counts)
}




# 建立空的結果表
max_summary <- data.frame(r_value = character(),
                          max_pathway = character(),
                          max_value = numeric(),
                          stringsAsFactors = FALSE)

# 針對每個 r_* 執行
for (r_name in names(pathway_counts_all)) {
  counts <- unlist(pathway_counts_all[[r_name]])
  max_val <- max(counts)
  max_paths <- names(counts)[counts == max_val]
  
  # 把每個 pathway 都記錄進結果表
  for (p in max_paths) {
    max_summary <- rbind(max_summary,
                         data.frame(r_value = r_name,
                                    max_pathway = p,
                                    max_value = max_val))
  }
}

print(max_summary)


high_pathways <- list()

for (r in r_values){
  r_name <- paste0("r_", r)
  count <- pathway_counts_all[[r_name]]
  high_pathways[[r_name]] <- names(count[count == 100])
}

intersection_pathways <- Reduce(intersect,high_pathways)

intersection_pathways

