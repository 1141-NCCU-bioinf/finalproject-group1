# ==============================================================================
# Project: RNA-seq Pathway Analysis - CTCF-KO Focus (V2: High Plot)
# ==============================================================================

# 1. 載入必要的套件
library(DESeq2)
library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)
library(stringr) # 用於文字換行

# 2. 讀取 Count Matrix
setwd("/Users/brian/Desktop/bio_final/rna") 

if (!file.exists("reconstructed_counts.csv")) {
  stop("錯誤: 找不到 reconstructed_counts.csv，請先執行 analysis.R！")
}

message("Loading count matrix...")
count_data <- read.csv("reconstructed_counts.csv", row.names = 1)

# 3. 設定 DESeq2 實驗設計
condition <- factor(c("Control", "Control", "Cp190", "Cp190", "CTCF", "CTCF"),
                    levels = c("Control", "Cp190", "CTCF"))

colData <- data.frame(row.names = colnames(count_data), condition = condition)

# 4. 執行 DESeq2
message("Running DESeq2...")
dds <- DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = colData,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ] # 過濾低表達基因
dds <- DESeq(dds)

# 5. 提取 CTCF-KO vs Control 的結果
message("Extracting results for CTCF-KO...")
res_CTCF <- results(dds, contrast=c("condition", "CTCF", "Control"))

# 篩選顯著基因 (FDR < 0.05, |Log2FC| > 1)
sig_genes_ctcf <- rownames(res_CTCF[which(res_CTCF$padj < 0.05 & abs(res_CTCF$log2FoldChange) > 1), ])
message("Significant genes found (CTCF vs Control): ", length(sig_genes_ctcf))

# 6. 執行 GO Enrichment 分析
if (length(sig_genes_ctcf) > 0) {
  message("Running GO Enrichment...")

  ego_ctcf <- enrichGO(gene          = sig_genes_ctcf,
                       OrgDb         = org.Dm.eg.db,
                       keyType       = "FLYBASE",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05)

  # 7. 畫圖並存檔 (改良版)
  if (!is.null(ego_ctcf) && nrow(ego_ctcf) > 0) {

    # 使用 stringr 來做文字換行 (width=50 表示每50個字元就換一行)
    p <- dotplot(ego_ctcf, showCategory=20) + 
      ggtitle("GO Enrichment: CTCF-KO vs Control") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) + 
      theme(axis.text.y = element_text(size = 10))

    # 設定高度 height = 12 (拉長圖片)
    ggsave("GO_Enrichment_CTCF_v2.png", p, width = 8, height = 12, dpi = 300)
    message("Success! Plot saved to GO_Enrichment_CTCF_v2.png")

  } else {
    message("No enriched GO terms found for CTCF-KO.")
  }
} else {
  message("No significant DE genes found for CTCF-KO.")
}
