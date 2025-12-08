# ==============================================================================
# 步驟 1: 安裝與載入套件
# ==============================================================================
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rtracklayer", "GenomicRanges", "GenomicFeatures", 
                       "TxDb.Dmelanogaster.UCSC.dm6.ensGene", "DESeq2", 
                       "clusterProfiler", "org.Dm.eg.db", "ggplot2"), update=FALSE, ask=FALSE)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(DESeq2)
library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)

# ==============================================================================
# 步驟 2: 準備基因註釋 (Gene Annotation)
# ==============================================================================
message("Loading dm6 gene annotation...")
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
exons_by_gene <- exonsBy(txdb, by="gene")

# ==============================================================================
# 步驟 3: 定義 Count 計算函數 (修正版)
# ==============================================================================
setwd("/Users/brian/Desktop/bio_final/rna") 
bed_files <- list.files(pattern = ".*dm6.bed.gz$")

if(length(bed_files) == 0) {
  stop("錯誤: 找不到 .bed.gz 檔案！請確認您的工作目錄。")
}

calc_counts_from_bed <- function(file_path, gene_regions) {
  message("Processing: ", file_path)

  # 讀取 BED 檔 (不需要 score)
  bg <- import(file_path, format = "BED")

  # 找出落在基因區域內的片段
  # ignore.strand=TRUE 確保我們能抓到所有 reads，不論正負鏈標記
  hits <- findOverlaps(bg, gene_regions, ignore.strand=TRUE)

  # 直接計算每個基因被多少個 BED reads 覆蓋
  counts <- countSubjectHits(hits)

  # 構建完整向量
  res <- counts
  names(res) <- names(gene_regions)

  return(res)
}

# ==============================================================================
# 步驟 4: 執行計算
# ==============================================================================
samples <- unique(gsub("_str[12]_dm6.bed.gz", "", bed_files))
message("Detected samples: ", paste(samples, collapse=", "))

count_matrix <- matrix(0, nrow=length(exons_by_gene), ncol=length(samples))
rownames(count_matrix) <- names(exons_by_gene)
colnames(count_matrix) <- samples

for (samp in samples) {
  f1 <- paste0(samp, "_str1_dm6.bed.gz")
  f2 <- paste0(samp, "_str2_dm6.bed.gz")

  counts_total <- rep(0, length(exons_by_gene))

  if (file.exists(f1)) counts_total <- counts_total + calc_counts_from_bed(f1, exons_by_gene)
  if (file.exists(f2)) counts_total <- counts_total + calc_counts_from_bed(f2, exons_by_gene)

  count_matrix[, samp] <- counts_total
}

# 存檔檢查
write.csv(count_matrix, "reconstructed_counts.csv")
message("Count matrix saved. Head of matrix:")
print(head(count_matrix))

# ==============================================================================
# 步驟 5: 進行 DESeq2 分析
# ==============================================================================
condition <- factor(c("Control", "Control", "Cp190", "Cp190", "CTCF", "CTCF"),
                    levels = c("Control", "Cp190", "CTCF"))

colData <- data.frame(row.names = colnames(count_matrix), condition = condition)

# 檢查是否全為 0
if (sum(count_matrix) == 0) {
  stop("錯誤: Count Matrix 總和為 0！計算邏輯仍有問題。")
}

dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData = colData,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# 提取結果: Cp190-KO vs Control
res_Cp190 <- results(dds, contrast=c("condition", "Cp190", "Control"))

# ==============================================================================
# 步驟 6: Pathway Analysis (GO Enrichment)
# ==============================================================================
sig_genes <- rownames(res_Cp190[which(res_Cp190$padj < 0.05 & abs(res_Cp190$log2FoldChange) > 1), ])
message("Significant genes found: ", length(sig_genes))

if (length(sig_genes) > 0) {
  ego <- enrichGO(gene          = sig_genes,
                  OrgDb         = org.Dm.eg.db,
                  keyType       = "FLYBASE",  
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05)

  if (!is.null(ego) && nrow(ego) > 0) {
    p <- dotplot(ego, showCategory=15) + ggtitle("GO Enrichment: Cp190-KO vs Control")
    ggsave("GO_Enrichment_Cp190.png", p)
    message("Plot saved to GO_Enrichment_Cp190.png")
  } else {
    message("No enriched GO terms found.")
  }
}
