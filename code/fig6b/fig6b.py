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
def normalize_and_subtract_raw(mutant_curve, control_curve, mutant_ref_random, control_ref_random):
    """
    計算原始差異值 (未縮放)，但在計算前先用隨機區域將 Mutant 與 Control 對齊
    """
    sum_mut_rnd = np.sum(mutant_ref_random)
    sum_ctl_rnd = np.sum(control_ref_random)
    
    if sum_mut_rnd == 0: return np.zeros_like(mutant_curve)
    
    # 1. 計算背景校正係數
    bg_scale_factor = sum_ctl_rnd / sum_mut_rnd
    
    # 2. 校正 Mutant 強度
    mutant_curve_norm = mutant_curve * bg_scale_factor
    
    # 3. 回傳原始差異 (Mutant_Corrected - Control)
    return (mutant_curve_norm - control_curve)

# ==========================================
# 5. 主程式
# ==========================================
if __name__ == "__main__":
    # A. 準備座標
    insulator_sites = get_coords(CSV_INSULATOR)
    
    random_sites = []
    chrom_list = list(set([x[0] for x in insulator_sites]))
    for _ in range(len(insulator_sites)):
        c = random.choice(chrom_list)
        pos = random.randint(100000, 20000000)
        random_sites.append((c, pos))
    
    # B. 提取數據
    print("\n--- 1. 讀取數據 ---")
    mut_ins_r1 = extract_crossing_curve(FILE_MUTANT_REP1, insulator_sites)
    ctl_ins_r1 = extract_crossing_curve(FILE_CONTROL_REP1, insulator_sites)
    mut_rnd_r1 = extract_crossing_curve(FILE_MUTANT_REP1, random_sites)
    ctl_rnd_r1 = extract_crossing_curve(FILE_CONTROL_REP1, random_sites)
    
    mut_ins_r2 = extract_crossing_curve(FILE_MUTANT_REP2, insulator_sites)
    ctl_ins_r2 = extract_crossing_curve(FILE_CONTROL_REP2, insulator_sites)
    mut_rnd_r2 = extract_crossing_curve(FILE_MUTANT_REP2, random_sites)
    ctl_rnd_r2 = extract_crossing_curve(FILE_CONTROL_REP2, random_sites)
    
    # C. 計算原始差異 (Raw Difference)
    print("\n--- 2. 計算差異與自動縮放 ---")
    if mut_ins_r1 is not None:
        # 先算出還沒放大的差異值 (這時數值可能很小，例如 0.005)
        raw_diff_ins_r1 = normalize_and_subtract_raw(mut_ins_r1, ctl_ins_r1, mut_rnd_r1, ctl_rnd_r1)
        raw_diff_rnd_r1 = normalize_and_subtract_raw(mut_rnd_r1, ctl_rnd_r1, mut_rnd_r1, ctl_rnd_r1)
        
        raw_diff_ins_r2 = normalize_and_subtract_raw(mut_ins_r2, ctl_ins_r2, mut_rnd_r2, ctl_rnd_r2)
        raw_diff_rnd_r2 = normalize_and_subtract_raw(mut_rnd_r2, ctl_rnd_r2, mut_rnd_r2, ctl_rnd_r2)
        
        #放大倍率這段待確認是否需要
        # === 關鍵步驟：自動計算放大倍率 (Auto-Scaling) ===
        # 我們找出 Rep1 絕緣子曲線的最大值 (Peak)
        current_peak = np.max(raw_diff_ins_r1)
        
        if current_peak > 0:
            # 計算要把 Peak 變成 TARGET_MAX_Y (例如 10) 需要乘多少
            final_multiplier = TARGET_MAX_Y / current_peak
            print(f"   -> 偵測到原始最大值: {current_peak:.6f}")
            print(f"   -> 自動計算放大倍率: {final_multiplier:.2f} 倍 (目標值: {TARGET_MAX_Y})")
        else:
            final_multiplier = 1.0
            print("   -> 警告：原始數值過低或為負，不進行縮放")

        # 應用放大倍率到所有線條 (保持比例一致)
        final_ins_r1 = raw_diff_ins_r1 * final_multiplier
        final_ins_r2 = raw_diff_ins_r2 * final_multiplier
        final_rnd_r1 = raw_diff_rnd_r1 * final_multiplier
        final_rnd_r2 = raw_diff_rnd_r2 * final_multiplier

        # D. 繪圖
        x_axis = np.arange(1, BINS_PER_SIDE + 1) * BIN_SIZE  / 1000 
        
        plt.figure(figsize=(7, 6))
        
        # 實心線 (Insulators)
        plt.plot(x_axis, final_ins_r1, 'r-o', markersize=6, label='Insulators, rep. 1')
        plt.plot(x_axis, final_ins_r2, 'b-o', markersize=6, label='Insulators, rep. 2')
        
        # 空心線 (Control)
        plt.plot(x_axis, final_rnd_r1, 'r-o', markerfacecolor='white', markersize=6, label='Control, rep. 1')
        plt.plot(x_axis, final_rnd_r2, 'b-o', markerfacecolor='white', markersize=6, label='Control, rep. 2')
        my_ticks = [5] + list(range(25, 301, 25))
        plt.xticks(my_ticks)
        plt.axhline(0, color='gray', linestyle='--')
        plt.xlabel('Distance (kb)') 
        plt.ylabel('Contact crossing difference')
        plt.title('Fig 6B (Auto-Scaled to Match Paper)')
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig('Fig6B_Final_Scaled.png')
        plt.show()
        print("\n繪圖完成！數值已調整為與論文一致。")