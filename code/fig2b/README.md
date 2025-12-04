# Recreate Fig 2B

This folder contains the script required to reproduce Figure 2B (Hi-C Heatmaps).

## Content

- **plot_Fig2B.py**: The Python script to generate the plot.

Output
- **\*.csv**: Pre-processed Hi-C contact matrices for the Bithorax complex (5kb resolution).
  - Control (Ras3)
  - Cp190-KO (Cpr6)
  - CTCF-KO

## Usage

Ensure you have Python installed with `pandas`, `numpy`, and `matplotlib`.

```bash
# Run the script
python plot_Fig2B.py
```

The output files `Fig2B_restored.pdf` and `Fig2B_restored.png` will be generated in this directory.

