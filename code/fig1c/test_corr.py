import h5py
import numpy as np

# 只需要比較這兩個檔案
f1_path = "GSM5956403_RAS3_rep2_05k_dm6.gcmap"
f2_path = "GSM5956407_CTCF_rep2_05k_dm6.gcmap"

print(f"Comparing content of:\n1. {f1_path}\n2. {f2_path}")

try:
    with h5py.File(f1_path, 'r') as f1, h5py.File(f2_path, 'r') as f2:
        # 讀取一個染色體作為代表
        d1 = f1['chr2L/40kb'][:]
        d2 = f2['chr2L/40kb'][:]
        
        # 計算簡單相關性
        corr = np.corrcoef(d1.flatten(), d2.flatten())[0, 1]
        print(f"\nPearson Correlation (chr2L): {corr:.4f}")
        
        if corr > 0.98:
            print("\n!!! 警告: 這兩個檔案內容幾乎完全相同！ !!!")
            print("您可能將同一個檔案複製並重命名了，或者下載時出錯。")
        else:
            print("\n檔案內容不同 (這是好消息)。")
            
except Exception as e:
    print(f"Error: {e}")

