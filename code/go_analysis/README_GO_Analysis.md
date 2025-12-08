# GO Enrichment Analysis for Cp190 and CTCF Knockout

This repository documents the complete workflow for performing Gene Ontology (GO) Enrichment Analysis on RNA-seq data from the study *"Topological screen identifies hundreds of Cp190- and CTCF-dependent Drosophila chromatin insulator elements"* (Kahn et al., 2023).

The goal was to investigate the biological functions affected by the loss of **Cp190** and **CTCF** insulator proteins in *Drosophila melanogaster*.

## 1. Prerequisites & Data Source

### Data Source
The raw RNA-seq data was obtained from the Gene Expression Omnibus (GEO) under accession **GSE198761**.
- **Genomic Build**: dm6
- **Data Format**: BED files (Single-end reads mapped to the genome)

| Sample Name | Condition | GEO Accession | File Name (Example) |
| :--- | :--- | :--- | :--- |
| **Ras3_rep1** | Control | GSM5956408 | `GSM5956408_Ras3_rep1_str1_dm6.bed.gz` |
| **Ras3_rep2** | Control | GSM5956409 | `GSM5956409_Ras3_rep2_str1_dm6.bed.gz` |
| **CPR6_rep1** | Cp190-KO | GSM5956410 | `GSM5956410_CPR6_rep1_str1_dm6.bed.gz` |
| **CPR6_rep2** | Cp190-KO | GSM5956411 | `GSM5956411_CPR6_rep2_str1_dm6.bed.gz` |
| **CTCF_rep1** | CTCF-KO | GSM5956412 | `GSM5956412_CTCF_rep1_str1_dm6.bed.gz` |
| **CTCF_rep2** | CTCF-KO | GSM5956413 | `GSM5956413_CTCF_rep2_str1_dm6.bed.gz` |

*Note: Each sample has two strand-specific files (`str1`, `str2`) which were merged during analysis.*

### Environment & Tools
The analysis was performed using **R** on macOS.

**Required R Packages (Bioconductor):**
```r
BiocManager::install(c("rtracklayer", "GenomicRanges", "GenomicFeatures", 
                       "TxDb.Dmelanogaster.UCSC.dm6.ensGene", "DESeq2", 
                       "clusterProfiler", "org.Dm.eg.db", "ggplot2", "stringr"))
```

## 2. Methodology: From BED to Pathway

Since the provided raw data were BED files (genomic coordinates of mapped reads) rather than a standard Count Matrix, a custom reconstruction step was required.

### Workflow Overview
1.  **Count Matrix Reconstruction**: 
    - Loaded `dm6` gene annotations using `TxDb.Dmelanogaster.UCSC.dm6.ensGene`.
    - Imported BED files using `rtracklayer`.
    - Counted overlaps between BED reads and gene regions using `findOverlaps`.
    - Merged strand-specific counts (`str1 + str2`) for each sample.
2.  **Differential Expression Analysis**:
    - Used **DESeq2** to normalize counts and calculate differential expression.
    - Comparisons: **Cp190-KO vs Control** and **CTCF-KO vs Control**.
3.  **Pathway Analysis**:
    - Used **clusterProfiler** for GO Enrichment Analysis (Biological Process).
    - Mapped genes using FlyBase IDs (`keyType = "FLYBASE"`).

## 3. Step-by-Step Reproduction Guide

### Step 1: Download Data
Download the `GSE198761_RAW.tar` from GEO and extract it into your working directory (e.g., `rna/`). Ensure all `.bed.gz` files are present.

### Step 2: Reconstruct Count Matrix
Run the `analysis.R` script to generate the count matrix from raw BED files.

```bash
Rscript analysis.R
```

**Key Code Snippet (Count Reconstruction):**
```r
# Load BED file and count overlaps with gene regions
bg <- import(file_path, format = "BED")
hits <- findOverlaps(bg, gene_regions, ignore.strand=TRUE)
counts <- countSubjectHits(hits)
```

**Output**: `reconstructed_counts.csv`

### Step 3: Run GO Enrichment Analysis
Once the count matrix is ready, run the specific analysis scripts for Cp190 and CTCF. These scripts handle DESeq2, filtering, and visualization.

**For Cp190-KO:**
```bash
Rscript GO_analysis_Cp190.R
```

**For CTCF-KO:**
```bash
Rscript GO_analysis_CTCF.R
```

## 4. Troubleshooting & Solutions

### Issue 1: "All samples have 0 counts"
*   **Problem**: Initial scripts resulted in a count matrix full of zeros.
*   **Cause**: The BED files lacked a `score` column (which is typical for bedGraph but not for raw alignment BEDs), causing weighted count calculations to fail.
*   **Solution**: Switched from weighted counting to `countSubjectHits()`, which simply counts the number of reads overlapping each gene region. This successfully reconstructed the raw read counts.

### Issue 2: Chromosome Naming Mismatch
*   **Problem**: Potential mismatch between UCSC (`chr2L`) and Ensembl (`2L`) chromosome naming conventions.
*   **Solution**: Verified that both the BED files and the `TxDb` package used UCSC style (`chr2L`), ensuring correct overlap detection.

### Issue 3: Plot Readability
*   **Problem**: GO term labels were too long and overlapped with the y-axis, and the plot was too short, causing points to crowd.
*   **Solution**: 
    - Used `stringr::str_wrap` to automatically wrap long labels.
    - Increased the `height` parameter in `ggsave` to 12 inches to provide ample vertical space.

## 5. Results
The analysis generated two visualization plots:
- `GO_Enrichment_Cp190_v2.png`: Shows that Cp190 loss significantly affects **ion transport** and **sensory perception** pathways.
- `GO_Enrichment_CTCF_v2.png`: Shows that CTCF loss impacts developmental and signaling pathways, with distinct differences from Cp190.
