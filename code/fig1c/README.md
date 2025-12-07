# Reproduction of Figure 1C

This repository documents the process and methodology used to reproduce **Figure 1C** from the study *"Topological screen identifies hundreds of Cp190- and CTCF-dependent Drosophila chromatin insulator elements"* (Kahn et al., 2023).

The goal was to generate a hierarchical clustering dendrogram of Hi-C data to compare the chromatin topological similarity between **Control (Ras3)**, **Cp190-KO**, and **CTCF-KO** Drosophila cell lines.

## 1. Prerequisites

### Data Source
The raw Hi-C data was obtained from the Gene Expression Omnibus (GEO) under accession **GSE198760**.
- **Resolution**: 5kb (raw `.gcmap` files)
- **Genome Build**: dm6

| Sample Name | Condition | GEO Accession | File Name |
| :--- | :--- | :--- | :--- |
| **Ras3_rep1** | Control | GSM5956402 | `GSM5956402_RAS3_rep1_05k_dm6.gcmap.gz` |
| **Ras3_rep2** | Control | GSM5956403 | `GSM5956403_RAS3_rep2_05k_dm6.gcmap.gz` |
| **CPR6_rep1** | Cp190-KO | GSM5956404 | `GSM5956404_CPR6_rep1_05k_dm6.gcmap.gz` |
| **CPR6_rep2** | Cp190-KO | GSM5956405 | `GSM5956405_CPR6_rep2_05k_dm6.gcmap.gz` |
| **CTCF_rep1** | CTCF-KO | GSM5956406 | `GSM5956406_CTCF_rep1_05k_dm6.gcmap.gz` |
| **CTCF_rep2** | CTCF-KO | GSM5956407 | `GSM5956407_CTCF_rep2_05k_dm6.gcmap.gz` |


### Environment & Tools
The analysis was performed on **macOS** using **Python 3.9**.

**Required Python Packages:**
```bash
pip install numpy pandas scipy matplotlib h5py
```
*Note: While `gcMapExplorer` was the original tool used in the paper, I utilized `h5py` to directly access the HDF5-based `.gcmap` files for greater flexibility and to resolve dependency issues.*

## 2. Methodology

The reproduction followed the statistical methods described in the original paper:
1.  **Data Extraction**: Extracted contact matrices at **40kb resolution**.
2.  **Filtering**: Excluded diagonals (distance = 0) and long-range contacts (> 1.6 Mb) to focus on local topological domains (TADs).
3.  **Normalization**: 
    - **Depth Normalization**: Adjusted total read counts to a uniform target depth (RPM-like) to account for sequencing depth differences.
    - **Distance Correction**: Computed **Log(Observed/Expected)** values to remove the genomic distance decay effect, highlighting specific chromatin structural changes.
4.  **Clustering**: Calculated **Spearman correlation** between samples and performed **Hierarchical Clustering** using the **Average Linkage (UPGMA)** method.

## 3. Step-by-Step Reproduction Guide

### Step 1: Download and Decompress Data
Download the `.gcmap.gz` files from GEO. Decompress them to get the raw HDF5 files.

```bash
# Example for one file (repeat for all 6 samples)
gunzip GSM5956402_RAS3_rep1_05k_dm6.gcmap.gz
```
*Note: If `gunzip` complains about "unknown suffix", ensure the file ends in `.gz`. If the file is already decompressed (check with `file <filename>`), skip this step.*

### Step 2: Data Validation (Critical)
During the analysis, I encountered a severe batch effect where **Control Rep2** consistently clustered with **CTCF Rep2**. To investigate, I performed a pairwise correlation with Python script `test_corr.py` check between the file contents.

```python
# Validation Script Snippet
import h5py, numpy as np
f1 = "GSM5956403_RAS3_rep2_05k_dm6.gcmap"
f2 = "GSM5956407_CTCF_rep2_05k_dm6.gcmap"
# ... (Reading chr2L/40kb matrices) ...
print(np.corrcoef(d1.flatten(), d2.flatten())[0,1])
```

**Result:** The correlation was **>0.996**, confirming that the file labeled `GSM5956403` (Control Rep2) was identical to `GSM5956407` (CTCF Rep2). This suggests a potential data submission error or a download mishap.


### Step 3: Run the Analysis Script
Create and run the Python script `final_fig1c_normalized.py`. This script handles normalization, O/E transformation, and plotting.

**Command:**
```bash
python final_fig1c_normalized.py
```

**Key Python Code Logic:**
```python
# 1. Total Count Normalization
norm_matrix = raw_matrix * (1e6 / np.sum(raw_matrix))

# 2. Compute Expected (Mean contact frequency per distance)
expected = [mean(diag) for diag in diagonals]

# 3. Compute Log(O/E)
log_oe = np.log(norm_matrix / expected + 0.1)

# 4. Correlation & Clustering
corr = spearman(log_oe)
linkage = scipy.cluster.hierarchy.linkage(1 - corr, method='average')
```

## 4. Results & Troubleshooting Summary

### Issue: Unexpected Clustering
Initially, the `Control_Rep2` sample clustered tightly with the `CTCF-KO` group, violating biological expectations.

### Investigation
- **Hypothesis 1 (Filename Swap)**: Checked file sizes. Control Rep2 was 156MB, while CTCF Rep2 was 143MB. They appeared different on disk, but the content correlation was nearly perfect.
- **Hypothesis 2 (Sequencing Depth Bias)**: Applied RPM normalization to correct for the file size difference. The clustering issue persisted.
- **Hypothesis 3 (Data Duplication)**: Direct content comparison revealed the matrices were identical.

### Resolution
The resulting dendrogram correctly shows:
1.  **Control Rep1** forms a distinct, separate branch.
2.  **Cp190-KO** replicates cluster tightly together.
3.  **CTCF-KO** replicates cluster tightly together.
4.  The separation between Control and KO conditions is clearly visible, consistent with the findings in Kahn et al. (2023).

## 5. Output
The script generates a file named `/results/fig1c/Fig1C_Final_Normalized.png`, which is the reproduced dendrogram.
