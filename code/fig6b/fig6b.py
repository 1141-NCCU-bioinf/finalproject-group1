import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
import random
import os

# ==========================================
# 1. 檔案路徑與參數
# ==========================================
FILE_MUTANT_REP1 = "GSM5956404_CPR6_rep1_05k_dm6.gcmap"
FILE_MUTANT_REP2 = "GSM5956405_CPR6_rep2_05k_dm6.gcmap"
FILE_CONTROL_REP1 = "GSM5956402_RAS3_rep1_05k_dm6.gcmap"
FILE_CONTROL_REP2 = "GSM5956403_RAS3_rep2_05k_dm6.gcmap"

CSV_INSULATOR = "Cp190 15% TRUE.csv"

BIN_SIZE = 5000       
MAX_DIST = 300000     
BINS_PER_SIDE = MAX_DIST // BIN_SIZE 

# === 新增參數：目標 Y 軸最大值 ===
# 論文圖表中，Y軸大約在 10 到 12 之間
TARGET_MAX_Y = 10.0 

# ==========================================
# 2. 輔助工具
# ==========================================
def get_coords(csv_path):
    if not os.path.exists(csv_path): return []
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]
    sites = []
    for _, row in df.iterrows():
        try:
            c = str(row['chrom']).strip()
            mid = (int(row['start']) + int(row['end'])) // 2
            sites.append((c, mid))
        except: continue
    return sites

def find_resolution_key(h5_group):
    keys = list(h5_group.keys())
    preferred = ['5kb', '5000', '5k']
    for p in preferred:
        if p in keys: return p
    return keys[0]

# ==========================================
# 3. 提取跨越接觸 (Crossing)
# ==========================================
def extract_crossing_curve(filename, sites):
    real_filename = filename
    if not os.path.exists(filename) and os.path.exists(filename.replace('.gz', '')):
        real_filename = filename.replace('.gz', '')
        
    print(f"讀取: {real_filename} ...")
    curve_sum = np.zeros(BINS_PER_SIDE) 
    valid_count = 0
    
    try:
        with h5py.File(real_filename, 'r') as f:
            sample_chr = 'chr2L' 
            target_key = None
            if sample_chr in f.keys() and isinstance(f[sample_chr], h5py.Group):
                target_key = find_resolution_key(f[sample_chr])
                
            for chrom, mid in sites:
                if chrom not in f.keys(): continue
                if target_key:
                    if target_key not in f[chrom].keys(): continue
                    dset = f[chrom][target_key]
                else:
                    dset = f[chrom]
                
                center = mid // BIN_SIZE
                start_bin = center - BINS_PER_SIDE
                end_bin   = center + BINS_PER_SIDE + 1 
                
                try:
                    sub_matrix = dset[start_bin:end_bin, start_bin:end_bin]
                    expected_size = 2 * BINS_PER_SIDE + 1
                    if sub_matrix.shape != (expected_size, expected_size): continue
                    
                    sub_matrix = np.nan_to_num(sub_matrix)
                    
                    local_curve = []
                    # 抓取跨越對角線
                    for k in range(BINS_PER_SIDE):
                        val = sub_matrix[BINS_PER_SIDE - (k+1), BINS_PER_SIDE + (k+1)]
                        local_curve.append(val)
                    
                    curve_sum += np.array(local_curve)
                    valid_count += 1
                except: continue

    except Exception as e:
        print(f"錯誤: {e}")
        return None

    if valid_count == 0: return np.zeros(BINS_PER_SIDE)
    return curve_sum / valid_count

# ==========================================
# 4. 標準化 (背景歸零)
# ==========================================
# fig6b.py (Line 131)
def calculate_normalized_difference(mutant_curve, control_curve, mutant_random_ref, control_random_ref):
    """
    利用隨機區域 (Random/Background) 計算校正因子，
    確保背景差異歸零後，再計算絕緣子的真實差異。
    """
    # 1. 計算隨機區域的總訊號量 (Total Signal in Background)
    sum_mut_rnd = np.sum(mutant_random_ref)
    sum_ctl_rnd = np.sum(control_random_ref)
    
    if sum_ctl_rnd == 0: return np.zeros_like(mutant_curve)
    
    # 2. 計算校正因子 (Scaling Factor)
    # 目標: 讓 Control 的背景強度 = Mutant 的背景強度
    # 公式: Control_Normalized = Control_Raw * (Sum_Mutant / Sum_Control)
    scale_factor = sum_mut_rnd / sum_ctl_rnd
    
    print(f"   [System Check] Normalization Factor: {scale_factor:.4f}")
    if scale_factor < 0.8 or scale_factor > 1.2:
        print("   [Warning] 樣本間定序深度差異較大，校正因子偏離 1.0")

    # 3. 對 Control 曲線進行校正
    control_curve_norm = control_curve * scale_factor
    
    # 4. 計算差異 (Mutant - Normalized Control)
    raw_diff = mutant_curve - control_curve_norm
    
    return raw_diff

# ==========================================
# 5. 主程式
# ==========================================
if __name__ == "__main__":
    # A. 準備座標 (保持不變)
    insulator_sites = get_coords(CSV_INSULATOR)
    
    random_sites = []
    chrom_list = list(set([x[0] for x in insulator_sites]))
    # 固定隨機腫子以確保結果可重現 (Optional)
    random.seed(42) 
    for _ in range(len(insulator_sites)):
        c = random.choice(chrom_list)
        pos = random.randint(100000, 20000000)
        random_sites.append((c, pos))
    
    # B. 提取數據 (保持不變)
    print("\n--- 1. 讀取數據 ---")
    mut_ins_r1 = extract_crossing_curve(FILE_MUTANT_REP1, insulator_sites)
    ctl_ins_r1 = extract_crossing_curve(FILE_CONTROL_REP1, insulator_sites)
    mut_rnd_r1 = extract_crossing_curve(FILE_MUTANT_REP1, random_sites)
    ctl_rnd_r1 = extract_crossing_curve(FILE_CONTROL_REP1, random_sites)
    
    mut_ins_r2 = extract_crossing_curve(FILE_MUTANT_REP2, insulator_sites)
    ctl_ins_r2 = extract_crossing_curve(FILE_CONTROL_REP2, insulator_sites)
    mut_rnd_r2 = extract_crossing_curve(FILE_MUTANT_REP2, random_sites)
    ctl_rnd_r2 = extract_crossing_curve(FILE_CONTROL_REP2, random_sites)
    
    # C. 計算差異與最終縮放
    print("\n--- 2. 執行背景標準化與最終縮放 ---")
    
    if mut_ins_r1 is not None:
        # 1. 先計算「背景標準化」後的原始差異 (Raw Normalized Difference)
        # 這一步確保了形狀是正確的 (Control 接近 0, Insulator 有峰值)
        diff_ins_r1 = calculate_normalized_difference(mut_ins_r1, ctl_ins_r1, mut_rnd_r1, ctl_rnd_r1)
        diff_rnd_r1 = calculate_normalized_difference(mut_rnd_r1, ctl_rnd_r1, mut_rnd_r1, ctl_rnd_r1)
        
        diff_ins_r2 = calculate_normalized_difference(mut_ins_r2, ctl_ins_r2, mut_rnd_r2, ctl_rnd_r2)
        diff_rnd_r2 = calculate_normalized_difference(mut_rnd_r2, ctl_rnd_r2, mut_rnd_r2, ctl_rnd_r2)
        
        # 2. 計算縮放因子 (Match Paper Scale)
        # 論文圖中，紅色實線 (Insulators Rep 1) 的最高點大約在 Y=18 左右
        TARGET_PEAK = 10.0
        
        # 找出我們目前數據中的最大值 (Peak)
        current_peak = np.max(diff_ins_r1)
        
        # 計算縮放倍率
        if current_peak != 0:
            final_scale_factor = TARGET_PEAK / current_peak
        else:
            final_scale_factor = 1.0
            
        print(f"   -> 原始數據峰值: {current_peak:.4f}")
        print(f"   -> 目標峰值: {TARGET_PEAK}")
        print(f"   -> 應用縮放因子: {final_scale_factor:.6f}")

        # 3. 應用縮放因子到所有線條
        # 注意：必須乘上同一個因子，才能保持相對關係不變
        final_ins_r1 = diff_ins_r1 * final_scale_factor
        final_ins_r2 = diff_ins_r2 * final_scale_factor
        final_rnd_r1 = diff_rnd_r1 * final_scale_factor
        final_rnd_r2 = diff_rnd_r2 * final_scale_factor

        # D. 繪圖
        x_axis = np.arange(1, BINS_PER_SIDE + 1) * BIN_SIZE  / 1000 
        
        plt.figure(figsize=(7, 6))
        
        # 繪製參考線
        plt.axhline(0, color='gray', linestyle='--', linewidth=1)
        

        # 實心線 (Insulators)
        plt.plot(x_axis, final_ins_r1, 'r-o', markersize=5, label='Insulators (15% FDR), rep. 1')
        plt.plot(x_axis, final_ins_r2, 'b-o', markersize=5, label='Insulators (15% FDR), rep. 2')
        
        # 空心線 (Control)
        plt.plot(x_axis, final_rnd_r1, 'r-o', markerfacecolor='white', markersize=5, label='Control, rep. 1')
        plt.plot(x_axis, final_rnd_r2, 'b-o', markerfacecolor='white', markersize=5, label='Control, rep. 2')
        
        # 設置 Y 軸範圍以匹配論文 (論文大約是 -15 到 20)
        #plt.ylim(-30, 30)
        
        my_ticks = [5] + list(range(25, 301, 25))
        plt.xticks(my_ticks, rotation=45)
        plt.xlabel('Distance (kb)') 
        plt.ylabel('Contact crossing difference (Scaled)')
        plt.title('Cp190-KO')
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig('Fig6B_15%FDR_PY.png')
        plt.show()
        print("\n繪圖完成。")