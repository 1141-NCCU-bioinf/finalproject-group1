# === 最終整合與視覺完美復現版 ===
library(sva)
library(limma)
library(dendextend) # 必須載入此套件來處理旋轉

# 1. 設定樣本與染色體
sample_bases <- c(
  "GSM5956402_RAS3_rep1_05k_dm6", "GSM5956403_RAS3_rep2_05k_dm6",
  "GSM5956404_CPR6_rep1_05k_dm6", "GSM5956405_CPR6_rep2_05k_dm6",
  "GSM5956406_CTCF_rep1_05k_dm6", "GSM5956407_CTCF_rep2_05k_dm6"
)
labels <- c("Control_r1", "Control_r2", "Cp190_r1", "Cp190_r2", "CTCF_r1", "CTCF_r2")

target_chroms <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")

# 2. 讀取並合併全基因組數據
cat("正在讀取全基因組數據...\n")
combined_data <- NULL

for (i in 1:length(sample_bases)) {
  # cat(paste0("  讀取樣本: ", labels[i], "... ")) 
  # (為了版面整潔，註解掉逐行印出，保留總體進度)
  sample_genome_vec <- c()
  
  for (chrom in target_chroms) {
    fname <- paste0(sample_bases[i], "_", chrom, "_40k_matrix.txt")
    
    if(file.exists(fname)) {
      mat <- as.matrix(read.table(fname))
      mat[mat < 2] <- NA 
      
      idx <- which(upper.tri(mat), arr.ind = TRUE)
      dist <- abs(idx[,1] - idx[,2])
      mask <- dist > 0 & dist <= 40
      
      mat_fixed <- as.matrix(read.table(fname))
      valid_vals_fixed <- mat_fixed[idx[mask,]]
      
      sample_genome_vec <- c(sample_genome_vec, valid_vals_fixed)
    }
  }
  
  if (is.null(combined_data)) {
    combined_data <- data.frame(sample_genome_vec)
  } else {
    if(nrow(combined_data) == length(sample_genome_vec)) {
      combined_data <- cbind(combined_data, sample_genome_vec)
    }
  }
}
colnames(combined_data) <- labels

# 3. 準備 ComBat
cat("數據讀取完成。開始進行 ComBat 批次修正...\n")
log_data <- log2(combined_data + 1)
batch <- c(1, 2, 1, 2, 1, 2)
mod <- model.matrix(~as.factor(c("Control", "Control", "Cp190", "Cp190", "CTCF", "CTCF")))
combat_data <- ComBat(dat = as.matrix(log_data), batch = batch, mod = mod, par.prior = TRUE)

# 4. 計算相關係數與聚類
cat("計算 Spearman 相關係數...\n")
cor_mat <- cor(combat_data, method = "spearman")
dist_mat <- as.dist(1 - cor_mat)
hclust_res <- hclust(dist_mat, method = "average")

# ==========================================
# === 5. 最終視覺調整 (Visual Tweaking) ===
# ==========================================

cat("正在繪製最終圖表 (Figure_1C_Final_Correction_v2.png)...\n")

# A. 準備旋轉 (確保 Control -> Cp190 -> CTCF)
dend <- as.dendrogram(hclust_res)
target_order_ids <- c("Control_r1", "Control_r2", "Cp190_r1", "Cp190_r2", "CTCF_r1", "CTCF_r2")
dend_rotated <- rotate(dend, target_order_ids)

# B. 更新標籤名稱 (Format Labels)
new_labels <- c(
  "Control (rep. 1)", "Control (rep. 2)", 
  "Cp190-KO (rep. 1)", "Cp190-KO (rep. 2)", 
  "CTCF-KO (rep. 1)", "CTCF-KO (rep. 2)"
)
labels(dend_rotated) <- new_labels

# C. 轉回 hclust 物件 (讓 hang 參數生效)
hclust_final <- as.hclust(dend_rotated)

# D. 設定畫布
# 【修改 1】將 width 從 500 增加到 550，給長標題更多橫向空間
png("Figure_1C_Final_Correction_v2.png", width = 550, height = 700, res = 100)

# E. 設定邊界
# 順序：下、左、上、右
# 【修改 2】將上邊界 (第3個數值) 從 5 增加到 6，給標題更多縱向空間
par(mar = c(7, 6, 6, 2))

# F. 計算 Y 軸刻度
min_h <- min(hclust_final$height)
max_h <- max(hclust_final$height)
yticks <- seq(from = min_h, to = max_h, length.out = 5) 
ytick_labels <- sprintf("%.3f", yticks) 

# G. 繪圖
plot(hclust_final,
     hang = 0.05,       # 短垂線
     axes = FALSE,      
     ann = FALSE,       
     main = "",
     cex = 0.9)         

# H. 手動繪製 Y 軸
axis(side = 2, at = yticks, labels = ytick_labels, las = 0, lwd = 1.5)         

# I. 加入標籤與標題
title(ylab = expression("Branch length " * (1 - rho)), line = 4, cex.lab = 1.2)
# 字體大小維持 1.3
title(main = "Reproduction of Hi-C experiment clustering", cex.main = 1.3)

dev.off()
cat("完成！請查看檔案: Figure_1C_Final_Correction_v2.png\n")


##################################################################
# ==========================================
# === 6. 加分題：佐證圖表 (修正版) ===
# ==========================================

# 安裝並載入必要的繪圖套件 (如果還沒裝過)
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# --- 準備工作：統一將名稱改為論文格式 ---
# 建立一個新名稱的向量，順序必須跟原始數據 (Control_r1...CTCF_r2) 一致
formal_names <- c(
  "Control (rep. 1)", "Control (rep. 2)", 
  "Cp190-KO (rep. 1)", "Cp190-KO (rep. 2)", 
  "CTCF-KO (rep. 1)", "CTCF-KO (rep. 2)"
)

# 替換相關係數矩陣的名稱
cor_mat_fixed <- cor_mat
rownames(cor_mat_fixed) <- formal_names
colnames(cor_mat_fixed) <- formal_names

# 替換 PCA 數據的名稱
combat_data_fixed <- combat_data
colnames(combat_data_fixed) <- formal_names


# --- 圖表 A: 標準樣式熱力圖 (Standard Heatmap) ---
png("Figure_S1_Heatmap_Standard.png", width = 800, height = 600, res = 100)

# 使用 pheatmap 繪製
# display_numbers = TRUE: 直接把相關係數寫在格子上 (最直觀)
# cluster_rows/cols = TRUE: 讓它自動把相似的聚在一起 (跟樹狀圖邏輯一樣)
pheatmap(cor_mat_fixed, 
         display_numbers = TRUE,   # 顯示數值
         number_format = "%.3f",   # 小數點後三位
         fontsize_number = 12,     # 數字大小
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), # 藍(低)-白-紅(高) 標準配色
         border_color = "grey60",  # 格子邊框
         main = "Spearman Correlation Matrix",
         angle_col = 45)           # 欄位名稱旋轉 45 度避免重疊

dev.off()
cat("圖表 A (標準熱力圖) 已輸出: Figure_S1_Heatmap_Standard.png\n")


# --- 圖表 B: PCA 主成分分析 (修正標籤與邊界) ---
png("Figure_S2_PCA_Fixed.png", width = 700, height = 600, res = 100)

# 1. 計算 PCA
pca_res <- prcomp(t(combat_data_fixed)) # 使用改名後的數據
pc1_var <- round(summary(pca_res)$importance[2, 1] * 100, 1)
pc2_var <- round(summary(pca_res)$importance[2, 2] * 100, 1)

# 2. 設定繪圖範圍 (解決字被切掉的問題)
# 找出 X 和 Y 軸的最大最小值，然後往外擴張 20%
x_limits <- range(pca_res$x[,1]) * 1.3
y_limits <- range(pca_res$x[,2]) * 1.3

# 3. 設定顏色與形狀
# Control=藍, Cp190=紅, CTCF=綠
group_colors <- c("blue", "blue", "red", "red", "green", "green") 
# Rep1=圓形(19), Rep2=三角形(17) (實心圖案比較清楚)
batch_shapes <- c(19, 17, 19, 17, 19, 17) 

# 4. 畫圖 (主體)
plot(pca_res$x[,1], pca_res$x[,2], 
     col = group_colors, 
     pch = batch_shapes,
     cex = 2.5, # 點大一點
     xlab = paste0("PC1 (", pc1_var, "%) - Biological Difference"),
     ylab = paste0("PC2 (", pc2_var, "%)"),
     main = "PCA: Biological Grouping vs. Batch",
     xlim = x_limits, # 套用擴大的範圍
     ylim = y_limits)

# 5. 加入標籤 (使用正式名稱)
# pos = 3 (上方), pos = 1 (下方) -> 自動錯開
text(pca_res$x[,1], pca_res$x[,2], 
     labels = rownames(pca_res$x), 
     pos = c(3, 1, 3, 1, 3, 1), # 讓標籤一上一下跳動，避免重疊
     cex = 0.9, font = 2)

# 6. 加入圖例
legend("bottomright", 
       legend = c("Control", "Cp190-KO", "CTCF-KO", "Rep. 1", "Rep. 2"),
       col = c("blue", "red", "green", "black", "black"),
       pch = c(15, 15, 15, 19, 17),
       pt.cex = 1.5,
       bg = "white")

dev.off()
cat("圖表 B (PCA) 已輸出: Figure_S2_PCA_Fixed.png\n")

