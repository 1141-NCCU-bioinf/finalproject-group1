    import h5py
    import numpy as np
    import glob
    import os

    # 目標解析度名稱
    TARGET_KEY_NAME = '40kb'

    # 設定要提取的染色體白名單 (避免提取到奇怪的碎片 contigs)
    # 根據果蠅基因組 (dm3/dm6)，主要染色體如下：
    VALID_CHROMS = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 
                    '2L', '2R', '3L', '3R', '4', 'X']

    def process_gcmap_whole_genome(filename):
        print(f"=== 正在處理: {filename} ===")
        
        try:
            with h5py.File(filename, 'r') as f:
                # 取得檔案內所有染色體名稱
                all_keys = list(f.keys())
                
                # 過濾出主要染色體
                target_chroms = [k for k in all_keys if k in VALID_CHROMS]
                
                if not target_chroms:
                    print(f"  [警告] 找不到標準染色體 (chr2L...chrX)。現有鍵值: {all_keys}")
                    # 如果沒有標準命名，則提取所有長度合理的鍵值作為備案
                    target_chroms = all_keys

                print(f"  -> 將提取以下染色體: {target_chroms}")
                
                extracted_count = 0
                
                for chrom in target_chroms:
                    # 檢查該染色體下是否有 40kb 數據
                    if TARGET_KEY_NAME in f[chrom]:
                        dataset_path = f"{chrom}/{TARGET_KEY_NAME}"
                        
                        # 讀取數據
                        data = f[dataset_path][:]
                        
                        # 建立輸出檔名：原檔名_染色體_40k_matrix.txt
                        # 例如: GSM5956402_..._chr2L_40k_matrix.txt
                        clean_fname = filename.replace('.gcmap', '')
                        output_filename = f"{clean_fname}_{chrom}_40k_matrix.txt"
                        
                        # 存檔
                        np.savetxt(output_filename, data, fmt='%g', delimiter='\t')
                        print(f"    [OK] {chrom} ({data.shape}) -> {output_filename}")
                        extracted_count += 1
                    else:
                        print(f"    [Skip] {chrom} : 找不到 '{TARGET_KEY_NAME}' 層")
                
                if extracted_count == 0:
                    print("  [失敗] 此檔案未提取出任何有效數據。")

        except Exception as e:
            print(f"  [錯誤] 處理檔案時發生例外: {e}")
        
        print("\n")

    # 主程式
    files = glob.glob("*.gcmap")
    if not files:
        print("找不到 .gcmap 檔案。")
    else:
        print(f"找到 {len(files)} 個 .gcmap 檔案，開始全基因組提取...\n")
        for f in files:
            process_gcmap_whole_genome(f)
        print("所有作業完成。")