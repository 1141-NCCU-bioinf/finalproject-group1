import h5py
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import os

# ================= 檔案清單 =================
FILES = {
    "Control_Rep1": "GSM5956402_RAS3_rep1_05k_dm6.gcmap",
    "Control_Rep2": "GSM5956403_RAS3_rep2_05k_dm6.gcmap",
    "Cp190_Rep1":   "GSM5956404_CPR6_rep1_05k_dm6.gcmap",
    "Cp190_Rep2":   "GSM5956405_CPR6_rep2_05k_dm6.gcmap",
    "CTCF_Rep1":    "GSM5956406_CTCF_rep1_05k_dm6.gcmap",
    "CTCF_Rep2":    "GSM5956407_CTCF_rep2_05k_dm6.gcmap"
}

CHROMS = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX']
RES_KEY = '40kb'
MAX_DIST_BINS = 40  # 1.6 Mb

# ================= 函數 =================

def compute_normalized_log_oe(matrix, max_bins, target_depth=1e6):
    """
    1. Total Count Normalization (RPM-like)
    2. Compute Log(O/E)
    """
    # 1. 深度標準化 (Depth Normalization)
    total_reads = np.nansum(matrix)
    if total_reads == 0: return None
    
    # 將矩陣標準化為 "每百萬 reads" (RPM) 或其他固定常數
    norm_matrix = matrix * (target_depth / total_reads)
    
    # 2. 計算 Expected (使用標準化後的矩陣)
    rows, cols = np.indices(norm_matrix.shape)
    dist_mat = np.abs(rows - cols)
    
    expected_by_dist = {}
    for d in range(1, max_bins + 1):
        mask = (dist_mat == d)
        vals = norm_matrix[mask]
        if len(vals) > 0:
            expected_by_dist[d] = np.mean(vals)
        else:
            expected_by_dist[d] = 1.0 # 避免除以0
            
    # 3. 計算 O/E
    mask_all = (dist_mat > 0) & (dist_mat <= max_bins)
    obs = norm_matrix[mask_all]
    dists = dist_mat[mask_all]
    exp = np.array([expected_by_dist.get(d, 1.0) for d in dists])
    
    # O/E
    oe = obs / (exp + 1e-9)
    
    # Log transform
    return np.log(oe + 0.1)

def load_data(filepath, chroms, res_key, max_bins):
    all_vals = []
    if not os.path.exists(filepath): return None
    
    try:
        with h5py.File(filepath, 'r') as f:
            for chrom in chroms:
                path = f"{chrom}/{res_key}"
                if path in f:
                    mat = f[path][:]
                    # 使用常數 target_depth=1e6 確保所有樣本縮放到同一水平
                    val = compute_normalized_log_oe(mat, max_bins, target_depth=1e6)
                    if val is not None:
                        all_vals.append(val)
    except Exception as e:
        print(f"Error {filepath}: {e}")
        return None
    
    if not all_vals: return None
    return np.concatenate(all_vals)

# ================= 執行 =================
print("Starting analysis with Depth Normalization...")
data = {}

for name, fname in FILES.items():
    print(f"Processing {name}...")
    vec = load_data(fname, CHROMS, RES_KEY, MAX_DIST_BINS)
    if vec is not None:
        data[name] = vec

if len(data) >= 3:
    # 確保長度一致
    min_len = min(len(v) for v in data.values())
    df = pd.DataFrame({k: v[:min_len] for k, v in data.items()})
    
    # 相關性與分群
    corr = df.corr(method='spearman')
    print("\nCorrelation Matrix:")
    print(corr)
    
    plt.figure(figsize=(10, 8))
    dist = 1 - corr
    linkage = sch.linkage(sch.distance.squareform(dist), method='average')
    
    sch.dendrogram(linkage, labels=corr.columns, leaf_rotation=45)
    plt.title("Reproduction of Fig 1C (Normalized Log O/E)")
    plt.ylabel("Distance (1 - Spearman)")
    plt.tight_layout()
    plt.savefig("Fig1C_Final_Normalized.png")
    plt.show()
    print("Done. Check Fig1C_Final_Normalized.png")
else:
    print("Not enough samples.")

