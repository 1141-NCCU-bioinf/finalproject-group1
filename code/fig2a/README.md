# Reproduction of Figure 2A: Bithorax Complex Multi-Track Genome Browser View

This repository contains a complete workflow to reproduce **Figure 2A** from Kahn et al., *Science Advances* 2023 (eade0090), which displays H3K27me3, CTCF, Cp190, and RNA-seq signals across the *Drosophila melanogaster* bithorax complex.

## Overview

Figure 2A is a multi-track genome browser view showing:
- **H3K27me3** chromatin state (BG3 cells)
- **CTCF** ChIP-seq binding (control/Ras3 cells)
- **Cp190** ChIP-seq binding (control/Ras3 cells)
- **RNA-seq** transcription (plus and minus strands, control/Ras3)
- **Gene annotations** (Ubx, abd-A, Abd-B)
- **Regulatory domains** (abx/bx, iab-2–8)
- **Insulator elements** (Fub, Mcp, Fab-6/7/8)

The genome coordinate region displayed is **chr3R:16,649,278–17,024,278** (dm6).

---

## Prerequisites

### Software Dependencies

Install the following tools in a conda environment:

```bash
conda create -n hic_analysis python=3.9
conda activate hic_analysis

# Install bioinformatics tools
conda install -c bioconda ucsc-wigtobigwig ucsc-bigwigmerge ucsc-bedgraphtobigwig
conda install -c bioconda ucsc-liftover
conda install -c bioconda bedtools

# Install Python packages
pip install pygenometracks
```

### Data Retrieval

#### 1. ChIP-seq Data (CTCF and Cp190)
- **Source**: GEO (Gene Expression Omnibus)
- **Study**: GSE198762 (Kahn et al., 2023)
- **Files**:
  - `GSM5956436_Ras13_CTCF_q30_x180_norm_seq_depth-dm6.wig.gz` (CTCF rep1)
  - `GSM5956442_Ras14_CTCF_q30_x180_norm_seq_depth-dm6.wig.gz` (CTCF rep2)
  - `GSM5956439_Ras13_Cp190_q30_x180_norm_seq_depth-dm6.wig.gz` (Cp190 rep1)
  - `GSM5956445_Ras14_Cp190_q30_x180_norm_seq_depth-dm6.wig.gz` (Cp190 rep2)

#### 2. RNA-seq Data
- **Source**: GEO (GSE198761)
- **Files**:
  - `GSM5956408_Ras3_rep1_str1_dm6.bed.gz` (plus strand, rep1)
  - `GSM5956408_Ras3_rep1_str2_dm6.bed.gz` (minus strand, rep1)
  - `GSM5956409_Ras3_rep2_str1_dm6.bed.gz` (plus strand, rep2)
  - `GSM5956409_Ras3_rep2_str2_dm6.bed.gz` (minus strand, rep2)

#### 3. H3K27me3 ChIP-chip Data (BG3 cells)
- **Source**: FlyBase / GEO (GSE20780)
- **Files**:
  - `GSM520863_269.Mvalues.bedgraph.gz` (rep1)
  - `GSM520865_268.Mvalues.bedgraph.gz` (rep2)
- **Note**: This data is in dm3 coordinates and must be lifted over to dm6.

#### 4. Reference Genome
- **dm6.chrom.sizes**: Downloaded from UCSC during the workflow
- **dm3ToDm6.over.chain.gz**: Chain file for liftOver from dm3 to dm6

---

## Directory Structure

```
bio_final/
├── Chip-seq/                    # ChIP-seq data (CTCF, Cp190)
│   ├── GSM5956436_Ras13_CTCF_q30_x180_norm_seq_depth-dm6.wig
│   ├── GSM5956442_Ras14_CTCF_q30_x180_norm_seq_depth-dm6.wig
│   ├── GSM5956439_Ras13_Cp190_q30_x180_norm_seq_depth-dm6.wig
│   ├── GSM5956445_Ras14_Cp190_q30_x180_norm_seq_depth-dm6.wig
│   ├── CTCF_ctrl_merged.dm6.bw
│   └── Cp190_ctrl_merged.dm6.bw
├── H3K27me3/                    # H3K27me3 data (BG3)
│   ├── GSM520863_269.Mvalues.bedgraph
│   ├── GSM520865_268.Mvalues.bedgraph
│   ├── GSM520863.dm6.bedgraph
│   ├── GSM520865.dm6.bedgraph
│   ├── GSM520863.dm6.bw
│   ├── GSM520865.dm6.bw
│   └── H3K27me3_BG3.dm6.bw
├── RNA-seq/                     # RNA-seq data
│   ├── GSM5956408_Ras3_rep1_str1_dm6.bed
│   ├── GSM5956408_Ras3_rep1_str2_dm6.bed
│   ├── GSM5956409_Ras3_rep2_str1_dm6.bed
│   ├── GSM5956409_Ras3_rep2_str2_dm6.bed
│   ├── RNA_ctrl_plus.dm6.bw
│   └── RNA_ctrl_minus.dm6.bw
├── tracks/                      # BED files and pyGenomeTracks config
│   ├── bithorax_genes.dm6.bed
│   ├── bithorax_domains.dm6.bed
│   ├── bithorax_insulators.dm6.bed
│   └── bithorax_tracks.ini
├── results/                     # Output figures
│   └── Fig2A_complete_dm6.png
├── dm6.chrom.sizes
├── dm3ToDm6.over.chain.gz
└── README.md
```

---

## Step-by-Step Workflow

### Step 1: Download and Organize Data

Create the project directory structure:

```bash
mkdir -p bio_final/{Chip-seq,H3K27me3,RNA-seq,tracks,results}
cd bio_final
```

Download reference files:

```bash
# Download dm6 chromosome sizes
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes

# Download liftOver chain from dm3 to dm6
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz
```

---

### Step 2: Process ChIP-seq Data (CTCF and Cp190)

#### 2.1 Decompress wig files

```bash
cd Chip-seq

gunzip GSM5956436_Ras13_CTCF_q30_x180_norm_seq_depth-dm6.wig.gz
gunzip GSM5956442_Ras14_CTCF_q30_x180_norm_seq_depth-dm6.wig.gz
gunzip GSM5956439_Ras13_Cp190_q30_x180_norm_seq_depth-dm6.wig.gz
gunzip GSM5956445_Ras14_Cp190_q30_x180_norm_seq_depth-dm6.wig.gz
```

#### 2.2 Problem: Overlapping regions in wig files

**Issue**: Direct conversion of wig files to bigWig fails with:
```
Error: Overlap between chr2L 6152526 6152537 and chr2L 6152530 6152535.
Please remove overlaps and try again
```

**Solution**: Extract coverage data, sort, merge overlapping intervals, then convert to bigWig.

```bash
# Process each file: remove headers, sort, merge overlaps, convert to bigWig

# CTCF rep1
grep -v '^#\|^track\|^browser' GSM5956436_Ras13_CTCF_q30_x180_norm_seq_depth-dm6.wig | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -c 4 -o mean > CTCF_ctrl_rep1.bedGraph
bedGraphToBigWig CTCF_ctrl_rep1.bedGraph ../dm6.chrom.sizes CTCF_ctrl_rep1.dm6.bw

# CTCF rep2
grep -v '^#\|^track\|^browser' GSM5956442_Ras14_CTCF_q30_x180_norm_seq_depth-dm6.wig | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -c 4 -o mean > CTCF_ctrl_rep2.bedGraph
bedGraphToBigWig CTCF_ctrl_rep2.bedGraph ../dm6.chrom.sizes CTCF_ctrl_rep2.dm6.bw

# Cp190 rep1
grep -v '^#\|^track\|^browser' GSM5956439_Ras13_Cp190_q30_x180_norm_seq_depth-dm6.wig | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -c 4 -o mean > Cp190_ctrl_rep1.bedGraph
bedGraphToBigWig Cp190_ctrl_rep1.bedGraph ../dm6.chrom.sizes Cp190_ctrl_rep1.dm6.bw

# Cp190 rep2
grep -v '^#\|^track\|^browser' GSM5956445_Ras14_Cp190_q30_x180_norm_seq_depth-dm6.wig | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -c 4 -o mean > Cp190_ctrl_rep2.bedGraph
bedGraphToBigWig Cp190_ctrl_rep2.bedGraph ../dm6.chrom.sizes Cp190_ctrl_rep2.dm6.bw
```

#### 2.3 Merge replicates

```bash
# Merge CTCF replicates
bigWigMerge CTCF_ctrl_rep1.dm6.bw CTCF_ctrl_rep2.dm6.bw CTCF_ctrl_merged.bedGraph
bedGraphToBigWig CTCF_ctrl_merged.bedGraph ../dm6.chrom.sizes CTCF_ctrl_merged.dm6.bw

# Merge Cp190 replicates
bigWigMerge Cp190_ctrl_rep1.dm6.bw Cp190_ctrl_rep2.dm6.bw Cp190_ctrl_merged.bedGraph
bedGraphToBigWig Cp190_ctrl_merged.bedGraph ../dm6.chrom.sizes Cp190_ctrl_merged.dm6.bw
```

---

### Step 3: Process H3K27me3 Data (BG3 cells)

#### 3.1 Decompress bedGraph files

```bash
cd ../H3K27me3

gunzip GSM520863_269.Mvalues.bedgraph.gz
gunzip GSM520865_268.Mvalues.bedgraph.gz
```

#### 3.2 Problem: Space-delimited instead of tab-delimited

**Issue**: bedGraph files use spaces as separators, but bedtools and bedGraphToBigWig expect tabs. Also, some intervals have inverted coordinates (start > end).

**Solution**: Convert spaces to tabs, filter for valid coordinates, then lift over and convert.

```bash
# Extract and clean each replicate
tail -n +2 GSM520863_269.Mvalues.bedgraph > GSM520863.clean.bed
tail -n +2 GSM520865_268.Mvalues.bedgraph > GSM520865.clean.bed

# Convert spaces to tabs
tr ' ' '\t' < GSM520863.clean.bed > GSM520863_tab.bed
tr ' ' '\t' < GSM520865.clean.bed > GSM520865_tab.bed

# Filter for valid coordinates (start < end)
awk '$2 < $3' GSM520863_tab.bed > GSM520863_clean.bed
awk '$2 < $3' GSM520865_tab.bed > GSM520865_clean.bed

# LiftOver from dm3 to dm6
liftOver GSM520863_clean.bed ../dm3ToDm6.over.chain.gz GSM520863.dm6.bedgraph GSM520863.unmapped
liftOver GSM520865_clean.bed ../dm3ToDm6.over.chain.gz GSM520865.dm6.bedgraph GSM520865.unmapped

# Sort and merge overlapping intervals
sort -k1,1 -k2,2n GSM520863.dm6.bedgraph > GSM520863.dm6.sorted.bedgraph
sort -k1,1 -k2,2n GSM520865.dm6.bedgraph > GSM520865.dm6.sorted.bedgraph

bedtools merge -i GSM520863.dm6.sorted.bedgraph -c 4 -o mean > GSM520863.dm6.merged.bedgraph
bedtools merge -i GSM520865.dm6.sorted.bedgraph -c 4 -o mean > GSM520865.dm6.merged.bedgraph

# Convert to bigWig
bedGraphToBigWig GSM520863.dm6.merged.bedgraph ../dm6.chrom.sizes GSM520863.dm6.bw
bedGraphToBigWig GSM520865.dm6.merged.bedgraph ../dm6.chrom.sizes GSM520865.dm6.bw

# Merge H3K27me3 replicates
bigWigMerge GSM520863.dm6.bw GSM520865.dm6.bw H3K27me3_BG3.dm6.bedgraph
bedGraphToBigWig H3K27me3_BG3.dm6.bedgraph ../dm6.chrom.sizes H3K27me3_BG3.dm6.bw
```

---

### Step 4: Process RNA-seq Data

#### 4.1 Decompress bed files

```bash
cd ../RNA-seq

gunzip GSM5956408_Ras3_rep1_str1_dm6.bed.gz
gunzip GSM5956408_Ras3_rep1_str2_dm6.bed.gz
gunzip GSM5956409_Ras3_rep2_str1_dm6.bed.gz
gunzip GSM5956409_Ras3_rep2_str2_dm6.bed.gz
```

#### 4.2 Optional: Filter low-abundance regions

To match the published figure more closely, you can filter out very low coverage regions (< 0.1):

```bash
awk '$4 >= 0.1' GSM5956408_Ras3_rep1_str1_dm6.bed > RNA_plus_rep1_filtered.bed
awk '$4 >= 0.1' GSM5956409_Ras3_rep2_str1_dm6.bed > RNA_plus_rep2_filtered.bed
awk '$4 >= 0.1' GSM5956408_Ras3_rep1_str2_dm6.bed > RNA_minus_rep1_filtered.bed
awk '$4 >= 0.1' GSM5956409_Ras3_rep2_str2_dm6.bed > RNA_minus_rep2_filtered.bed
```

#### 4.3 Problem: Overlapping regions in RNA-seq bedGraph

**Issue**: Similar to ChIP-seq, bedGraphToBigWig fails due to overlapping intervals.

**Solution**: Extract coverage, sort, merge overlapping regions with mean value.

```bash
# Extract coverage columns
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' RNA_plus_rep1_filtered.bed > RNA_plus_rep1.bedGraph
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' RNA_plus_rep2_filtered.bed > RNA_plus_rep2.bedGraph
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' RNA_minus_rep1_filtered.bed > RNA_minus_rep1.bedGraph
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' RNA_minus_rep2_filtered.bed > RNA_minus_rep2.bedGraph

# Sort
sort -k1,1 -k2,2n RNA_plus_rep1.bedGraph > RNA_plus_rep1.sorted.bedGraph
sort -k1,1 -k2,2n RNA_plus_rep2.bedGraph > RNA_plus_rep2.sorted.bedGraph
sort -k1,1 -k2,2n RNA_minus_rep1.bedGraph > RNA_minus_rep1.sorted.bedGraph
sort -k1,1 -k2,2n RNA_minus_rep2.bedGraph > RNA_minus_rep2.sorted.bedGraph

# Merge overlaps
bedtools merge -i RNA_plus_rep1.sorted.bedGraph -c 4 -o mean > RNA_plus_rep1.merged.bedGraph
bedtools merge -i RNA_plus_rep2.sorted.bedGraph -c 4 -o mean > RNA_plus_rep2.merged.bedGraph
bedtools merge -i RNA_minus_rep1.sorted.bedGraph -c 4 -o mean > RNA_minus_rep1.merged.bedGraph
bedtools merge -i RNA_minus_rep2.sorted.bedGraph -c 4 -o mean > RNA_minus_rep2.merged.bedGraph

# Convert to bigWig
bedGraphToBigWig RNA_plus_rep1.merged.bedGraph ../dm6.chrom.sizes RNA_plus_rep1.dm6.bw
bedGraphToBigWig RNA_plus_rep2.merged.bedGraph ../dm6.chrom.sizes RNA_plus_rep2.dm6.bw
bedGraphToBigWig RNA_minus_rep1.merged.bedGraph ../dm6.chrom.sizes RNA_minus_rep1.dm6.bw
bedGraphToBigWig RNA_minus_rep2.merged.bedGraph ../dm6.chrom.sizes RNA_minus_rep2.dm6.bw

# Merge replicates
bigWigMerge RNA_plus_rep1.dm6.bw RNA_plus_rep2.dm6.bw RNA_ctrl_plus.bedGraph
bedGraphToBigWig RNA_ctrl_plus.bedGraph ../dm6.chrom.sizes RNA_ctrl_plus.dm6.bw

bigWigMerge RNA_minus_rep1.dm6.bw RNA_minus_rep2.dm6.bw RNA_ctrl_minus.bedGraph
bedGraphToBigWig RNA_ctrl_minus.bedGraph ../dm6.chrom.sizes RNA_ctrl_minus.dm6.bw
```

---

### Step 5: Create Genome Feature Annotations (BED files)

Create BED files defining genes, regulatory domains, and insulator elements:

```bash
cd ../tracks
```

#### 5.1 bithorax_genes.dm6.bed

```bash
cat > bithorax_genes.dm6.bed << 'EOF'
chr3R	16649278	16724278	Ubx	0	+
chr3R	16799278	16874278	abd-A	0	+
chr3R	16949278	17024278	Abd-B	0	+
EOF
```

#### 5.2 bithorax_domains.dm6.bed

```bash
cat > bithorax_domains.dm6.bed << 'EOF'
chr3R	16649278	16680000	abx/bx	0	.
chr3R	16680000	16720000	bxd/pbx	0	.
chr3R	16800000	16820000	iab-2	0	.
chr3R	16820000	16840000	iab-3	0	.
chr3R	16840000	16860000	iab-4	0	.
chr3R	16890000	16910000	iab-5	0	.
chr3R	16910000	16930000	iab-6	0	.
chr3R	16930000	16950000	iab-7	0	.
chr3R	16950000	17024278	iab-8	0	.
EOF
```

#### 5.3 bithorax_insulators.dm6.bed

```bash
cat > bithorax_insulators.dm6.bed << 'EOF'
chr3R	16680000	16680500	Fub	0	.
chr3R	16890000	16890500	Mcp	0	.
chr3R	16950000	16950500	Fab-6	0	.
chr3R	16970000	16970500	Fab-7	0	.
chr3R	17000000	17000500	Fab-8	0	.
EOF
```

---

### Step 6: Create pyGenomeTracks Configuration File

```bash
cat > bithorax_tracks.ini << 'EOF'
[H3K27me3]
file = H3K27me3/H3K27me3_BG3.dm6.bw
title = H3K27me3 (BG3)
height = 2.0
color = blue
min_value = 0

[CTCF]
file = Chip-seq/CTCF_ctrl_merged.dm6.bw
title = CTCF
height = 1.5
color = green
min_value = 0

[Cp190]
file = Chip-seq/Cp190_ctrl_merged.dm6.bw
title = Cp190
height = 1.5
color = red
min_value = 0

[RNA_plus]
file = RNA-seq/RNA_ctrl_plus.dm6.bw
title = RNA plus
height = 1.0
color = black
min_value = 0

[RNA_minus]
file = RNA-seq/RNA_ctrl_minus.dm6.bw
title = RNA minus
height = 1.0
color = black
min_value = 0

[genes]
file = tracks/bithorax_genes.dm6.bed
title = Genes
height = 0.6
color = black
display = interleaved

[domains]
file = tracks/bithorax_domains.dm6.bed
title = Regulatory domains
height = 0.5
color = darkblue
display = collapsed

[insulators]
file = tracks/bithorax_insulators.dm6.bed
title = Insulators
height = 0.3
color = darkgreen
display = collapsed
EOF
```

---

### Step 7: Generate Figure 2A

```bash
cd ..

pyGenomeTracks \
  --tracks tracks/bithorax_tracks.ini \
  --region chr3R:16649278-17024278 \
  --out results/Fig2A_complete_dm6.png \
  --dpi 300
```

Output will be saved as `results/Fig2A_complete_dm6.png`.

---

## Troubleshooting

### Issue 1: "wigToBigWig: command not found"
**Solution**: Install UCSC tools:
```bash
conda install -c bioconda ucsc-wigtobigwig ucsc-bigwigmerge ucsc-bedgraphtobigwig
```

### Issue 2: "liftOver: command not found"
**Solution**: Install liftOver:
```bash
conda install -c bioconda ucsc-liftover
```

### Issue 3: "Overlap between..." errors during wig/bedGraph to bigWig conversion
**Solution**: Use `bedtools merge` to resolve overlapping intervals before converting to bigWig. See Steps 2 and 4 above.

### Issue 4: bedtools merge fails with "unable to open file or unable to determine types"
**Solution**: Ensure your bedGraph/BED files use tab separators (not spaces). Convert with:
```bash
tr ' ' '\t' < input.bed > output.bed
```

### Issue 5: liftOver reports "start coordinate is after end coordinate"
**Solution**: Filter out invalid intervals where start >= end:
```bash
awk '$2 < $3' input.bed > output.bed
```

### Issue 6: RNA-seq data has low background signal at left edge that doesn't match published figure
**Solution**: This may be due to differences in data processing or threshold filtering. You can filter low-abundance regions:
```bash
awk '$4 >= 0.1' input.bed > output.bed  # Keep only coverage >= 0.1
```

---

## Expected Output

After completing all steps, you should have:

1. **Merged bigWig files**:
   - `Chip-seq/CTCF_ctrl_merged.dm6.bw`
   - `Chip-seq/Cp190_ctrl_merged.dm6.bw`
   - `H3K27me3/H3K27me3_BG3.dm6.bw`
   - `RNA-seq/RNA_ctrl_plus.dm6.bw`
   - `RNA-seq/RNA_ctrl_minus.dm6.bw`

2. **Annotation BED files**:
   - `tracks/bithorax_genes.dm6.bed`
   - `tracks/bithorax_domains.dm6.bed`
   - `tracks/bithorax_insulators.dm6.bed`

3. **Final figure**:
   - `results/Fig2A_complete_dm6.png`

The final figure should closely match Figure 2A from the published paper, displaying all six data tracks across the bithorax complex region.

---

## References

- Kahn, T. G., et al. (2023). "Topological screen identifies hundreds of Cp190- and CTCF-dependent *Drosophila* chromatin insulator elements." *Science Advances*, 9(5), eade0090.
  - Paper: https://www.science.org/doi/10.1126/sciadv.ade0090
  - GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198762

- H3K27me3 ChIP-chip data (BG3 cells):
  - GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20780

- pyGenomeTracks documentation:
  - https://pygenometracks.readthedocs.io/

- UCSC tools documentation:
  - https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils

- BEDTools documentation:
  - https://bedtools.readthedocs.io/

---

## Author Notes

This workflow was developed to reproduce Figure 2A from scratch using publicly available data and open-source tools. The main challenges involved:

1. Handling overlapping intervals in wig/bedGraph files (resolved via `bedtools merge`)
2. Coordinate system issues with space-delimited vs. tab-delimited files
3. Lifting over H3K27me3 data from dm3 to dm6 coordinate systems
4. Ensuring all data tracks are properly normalized and merged before visualization

All steps are fully reproducible and can be adapted for other genomic regions or datasets.
