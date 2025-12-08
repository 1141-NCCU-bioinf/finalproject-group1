library(rtracklayer)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# 設定路徑
setwd("/Users/brian/Desktop/bio_final/rna")

# 讀取一個 BED 檔
bed_file <- list.files(pattern = ".*dm6.bed.gz$")[1]
message("Checking file: ", bed_file)

bg <- import(bed_file, format="BED")

# 1. 檢查 BED 的染色體名稱
message("BED Chromosomes (前 5 個):")
print(head(seqlevels(bg)))

# 2. 檢查 TxDb 的染色體名稱
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
message("TxDb Chromosomes (前 5 個):")
print(head(seqlevels(txdb)))

# 3. 檢查是否有 Score 欄位
message("BED Data Columns:")
print(colnames(mcols(bg)))

if (!is.null(mcols(bg)$score)) {
    message("Score Summary:")
    print(summary(mcols(bg)$score))
} else {
    message("No 'score' column found!")
}

# 4. 測試重疊
exons <- exonsBy(txdb, by="gene")
hits <- findOverlaps(bg, exons)
message("Total Overlaps found: ", length(hits))

