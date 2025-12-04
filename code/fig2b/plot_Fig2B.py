import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sys
import os

def plot_fig2b():
    print("Generating publication-quality Fig 2B (Deep Orange)...")
    
    files = {
        "Control": ["RAS3_dm6_bithorax_5kb.csv", "RAS3_rep2_dm6_bithorax_5kb.csv"],
        "Cp190-KO": ["CPR6_dm6_bithorax_5kb.csv", "CPR6_rep2_dm6_bithorax_5kb.csv"],
        "CTCF-KO": ["CTCF_dm6_bithorax_5kb.csv", "CTCF_rep2_dm6_bithorax_5kb.csv"]
    }
    
    # 1. Load Data & Sum
    matrices = {}
    coords = None
    
    for name, file_list in files.items():
        mat_sum = None
        for f in file_list:
            df = pd.read_csv(f, index_col=0)
            if coords is None: coords = df.index.astype(float).values
            mat = df.values.astype(float)
            if mat_sum is None: mat_sum = mat
            else: mat_sum += mat
        matrices[name] = mat_sum

    # 2. Color Scale (0 - 160) with Midpoint 40 and Brown End
    z_min = 0
    z_max = 160
    mid_point = 40 # Reverted to 40 for smoother gradient
    
    norm_mid = (mid_point - z_min) / (z_max - z_min)
    nodes = [0.0, norm_mid, 1.0]
    colors = ["white", "orange", "brown"] # Reverted to brown (Deep Orange)
    
    cmap = LinearSegmentedColormap.from_list("custom_hic", list(zip(nodes, colors)), N=256)
    
    # 3. Setup Plot
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), constrained_layout=True)
    
    # 4. Alignment
    tick_indices = np.arange(15, 91, 15)
    tick_labels = ["16649", "16724", "16799", "16874", "16949", "17024"]
    
    # Reference Lines
    bin_size = 5000
    start_bp = coords[0]
    ref_bp = [16799278, 16874278]
    ref_indices = [(bp - start_bp) / bin_size for bp in ref_bp]

    # 6. Plotting
    images = []
    for i, (name, mat) in enumerate(matrices.items()):
        ax = axes[i]
        img = ax.imshow(mat, cmap=cmap, vmin=z_min, vmax=z_max, 
                        interpolation='nearest', origin='lower')
        images.append(img)
        
        ax.set_title(name, fontsize=14, pad=10, fontweight='bold')
        
        ax.set_xticks(tick_indices)
        ax.set_xticklabels(tick_labels, fontsize=9)
        ax.set_xlabel("kb", fontsize=10)
        
        if i == 0:
            ax.set_yticks(tick_indices)
            ax.set_yticklabels(tick_labels, fontsize=9)
            ax.set_ylabel("kb", fontsize=10)
        else:
            ax.set_yticks([])
            
        for idx in ref_indices:
            ax.axvline(x=idx, color='#32CD32', linestyle='--', linewidth=1.5, alpha=0.9)
        
        for spine in ax.spines.values():
            spine.set_visible(True); spine.set_color('black'); spine.set_linewidth(1)

    # 7. Colorbar
    cbar = fig.colorbar(images[0], ax=axes, orientation='vertical', fraction=0.02, pad=0.02)
    cbar.set_label('Contact numbers', fontsize=10)
    cbar.set_ticks(np.arange(0, 161, 20))
    
    plt.savefig("Fig2B_restored.pdf", dpi=300, bbox_inches='tight')
    plt.savefig("Fig2B_restored.png", dpi=300, bbox_inches='tight')
    print("Done.")

if __name__ == "__main__":
    plot_fig2b()
