# Reproduction of Figure 2B: Hi-C Heatmaps

This directory contains the complete pipeline to reproduce **Figure 2B** from the raw Hi-C data (GSE198760).

The pipeline consists of two main steps:
1. **Data Extraction**: Extracting the specific genomic region (Bithorax Complex) from the raw genome-wide `.gcmap` files.
2. **Visualization**: Plotting the contact matrices as heatmaps using Python.

## 1. Prerequisites

You will need Python 3.10+ and the following libraries:

```bash
pip install pandas numpy matplotlib h5py
```

## 2. Data Preparation (Raw to CSV)

The raw Hi-C data is provided in `.gcmap` (HDF5) format. We use the script `extract_bithorax_gcmap.py` to extract the 5kb-resolution contact matrix for the Bithorax Complex region (`chr3R:16,649,278-17,024,278`).

The raw data is expected to be in `../GSE198760_RAW/`.

**Script:** `extract_bithorax_gcmap.py`

### Usage Example

To extract the matrix for one sample (e.g., Ras3 Control Rep 1):

```bash
python extract_bithorax_gcmap.py \
  --gcmap ../GSE198760_RAW/GSM5956402_Ras3_rep1_05k_dm6.gcmap.gz \
  --chrom 3R \
  --start 16649278 \
  --end 17024278 \
  --out RAS3_dm6_bithorax_5kb.csv
```

*Note: This step has already been performed for you. The resulting `.csv` files are included in results/fig2b.*

## 3. Visualization (CSV to Figure)

Once the `.csv` matrices are ready, we use `plot_Fig2B.py` to sum the replicates and generate the publication-quality figure.

**Script:** `plot_Fig2B.py`

### Usage

```bash
python plot_Fig2B.py
```

### Output
The script will generate:
- `Fig2B_restored.png` (Preview image)


