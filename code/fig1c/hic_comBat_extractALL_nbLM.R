# === 最終整合版：全基因組數據 + ComBat 批次修正 ===
library(sva)
library(limma)

# 1. 設定樣本與染色體
sample_bases <- c(
  "GSM5956402_RAS3_rep1_05k_dm6", "GSM5956403_RAS3_rep2_05k_dm6",
  "GSM5956404_CPR6_rep1_05k_dm6", "GSM5956405_CPR6_rep2_05k_dm6",
  "GSM5956406_CTCF_rep1_05k_dm6", "GSM5956407_CTCF_rep2_05k_dm6"
)
labels <- c("Control_r1", "Control_r2", "Cp190_r1", "Cp190_r2", "CTCF_r1", "CTCF_r2")

# 請確保您的 Python 腳本已經跑完，並且有產生這些染色體的檔案
target_chroms <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")

# 2. 讀取並合併全基因組數據
cat("正在讀取全基因組數據...\n")
combined_data <- NULL

for (i in 1:length(sample_bases)) {
  cat(paste0("  讀取樣本: ", labels[i], "... "))
  sample_genome_vec <- c()
  
  for (chrom in target_chroms) {
    # 嘗試讀取檔案
    fname <- paste0(sample_bases[i], "_", chrom, "_40k_matrix.txt")
    
    # 容錯處理 (如果有些染色體沒提取到)
    if(file.exists(fname)) {
      mat <- as.matrix(read.table(fname))
      
      # 論文過濾步驟
      mat[mat < 2] <- NA # 過濾 < 2 (設為 NA)
      
      # 取上三角且距離 <= 40 bins
      idx <- which(upper.tri(mat), arr.ind = TRUE)
      dist <- abs(idx[,1] - idx[,2])
      mask <- dist > 0 & dist <= 40
      
      valid_vals <- mat[idx[mask,]]
      
      # 移除 NA (ComBat 不喜歡 NA)
      # 為了保持數據完整性，我們這裡簡單移除 NA 點
      # 注意：這樣不同樣本長度會不同，無法合併。
      # 修正策略：我們必須讓所有樣本長度一致才能跑 ComBat。
      # 最好的方法是：不過濾成 NA，而是補 0 或補最小數值。
      # 為了 ComBat，我們將 <2 的值視為 1 (Log2 後為 0)
      
      mat_fixed <- as.matrix(read.table(fname))
      valid_vals_fixed <- mat_fixed[idx[mask,]]
      
      sample_genome_vec <- c(sample_genome_vec, valid_vals_fixed)
    }
  }
  
  cat(paste("數據點數量:", length(sample_genome_vec), "\n"))
  
  if (is.null(combined_data)) {
    combined_data <- data.frame(sample_genome_vec)
  } else {
    # 確保長度一致 (理論上 Python 提取的矩陣大小都一樣，所以這裡應該安全)
    if(nrow(combined_data) == length(sample_genome_vec)) {
      combined_data <- cbind(combined_data, sample_genome_vec)
    } else {
      cat("    [警告] 數據長度不一致，跳過此樣本。\n")
    }
  }
}
colnames(combined_data) <- labels

# 3. 準備 ComBat
cat("數據讀取完成。開始進行 ComBat 批次修正...\n")

# Log 轉換 (ComBat 需要常態分佈)
# 這裡我們加 1，這樣原本 <2 的數值就會變成很小的數，不會報錯
log_data <- log2(combined_data + 1)

# 設定批次 (Batch) 和 分組 (Model)
# Batch: 1=Rep1, 2=Rep2
batch <- c(1, 2, 1, 2, 1, 2)
# Group: Control, Cp190, CTCF
mod <- model.matrix(~as.factor(c("Control", "Control", "Cp190", "Cp190", "CTCF", "CTCF")))

# 執行 ComBat
combat_data <- ComBat(dat = as.matrix(log_data), batch = batch, mod = mod, par.prior = TRUE)

# 4. 繪圖 (使用 Spearman)
cat("計算 Spearman 相關係數...\n")
cor_mat <- cor(combat_data, method = "spearman")

# 輸出矩陣檢查 (看看 Control_r1 vs Control_r2 有沒有變高)
print(round(cor_mat, 3))

# 繪製樹狀圖
dist_mat <- as.dist(1 - cor_mat)
hclust_res <- hclust(dist_mat, method = "average")

# 5. (選擇性) 視覺旋轉，讓 Control 在左邊
dend <- as.dendrogram(hclust_res)
# 嘗試讓 Control 排在最前面
target_order <- c("Control_r1", "Control_r2", "Cp190_r1", "Cp190_r2", "CTCF_r1", "CTCF_r2")

# 檢查是否安裝 dendextend
if (require("dendextend")) {
  dend_rotated <- rotate(dend, target_order)
  plot(dend_rotated, 
       main = "Reproduction of Figure 1C (Whole Genome + ComBat)", 
       sub = "The Final Solution: Data Quantity + Batch Correction")
} else {
  plot(hclust_res, 
       main = "Reproduction of Figure 1C (Whole Genome + ComBat)", 
       sub = "The Final Solution: Data Quantity + Batch Correction")
}

