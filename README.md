# Gut Microbiota and Coronary Artery Calcification: MR Meta-Analysis

[![R](https://img.shields.io/badge/R-%3E%3D4.3-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reproducible R pipeline for the two-sample Mendelian Randomization meta-analysis:

> **Gut Microbiota Composition and the Risk of Coronary Artery Calcification:
> A Genome-Wide Mendelian Randomization Meta-Analysis**

---

## Overview

This pipeline identifies gut bacterial genera causally associated with coronary
artery calcification (CAC) using genome-wide two-sample MR, combining:

- **Exposure:** MiBioGen (N = 18,340; 16S rRNA; 131 genera) + Dutch Microbiome Project
  (N = 7,738; shotgun metagenomics; 207 taxa)
- **Outcome:** Meta-analysed CAC GWAS (GCST90278456 + GCST90503074; median
  effective N = 49,309)

### Key findings
- 81 genus-level taxa with Bonferroni-significant causal effects on CAC
- Lachnospiraceae + Ruminococcaceae systematically protective
  (mean OR = 0.988, p = 0.044; Mann-Whitney p = 0.011 vs other Firmicutes)
- Net mPGS R² = 0.627%; D10 vs D1 OR ≈ 1.55

---

## Directory structure

```
gut_microbiota_CAC_MR/
├── R/
│   ├── config.R                        # paths, thresholds (edit before running)
│   ├── 00_install_packages.R           # install all dependencies
│   ├── 01_prepare_data.R               # load instruments; meta-analyse CAC GWAS
│   ├── 02_run_mr.R                     # IVW, MR-Egger, WM; taxonomy audit
│   ├── 03_sensitivity_analyses.R       # Steiger, MVMR, LOO, funnel plots
│   ├── 04_taxonomic_pattern_analysis.R # family/group-level statistical tests
│   ├── 05_mpgs_prediction.R            # mPGS, R², AUC, PAF
│   └── 06_forest_plots.R              # publication-quality forest plots
├── data/
│   └── raw/                            # place input files here (see below)
├── results/                            # auto-created; all CSV/RDS outputs
├── figures/                            # auto-created; all PNG figures
├── run_pipeline.R                      # master script
└── README.md
```

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/Sete4321/gut_microbiota_CAC_MR.git
cd gut_microbiota_CAC_MR
```

### 2. Install R packages

```r
source("R/00_install_packages.R")
```

### 3. Download input data

Place the following files in `data/raw/`:

| File | Source | Access |
|------|--------|--------|
| `MBG_allHits_p1e4.txt` | MiBioGen GWAS (Kurilshikov et al., Nat Genet 2021) | IEU OpenGWAS: https://gwas.mrcieu.ac.uk or contact alexa.kur@gmail.com |
| `dmp_summary_stats_taxa.csv` | Dutch Microbiome Project (Lopera-Maya et al., Nature Genet 2022) | Dutch Microbiome Project: https://dutchmicrobiomeproject.molgeniscloud.org/ or contact serena.sanna@irgb.cnr.it |
| `GCST90278456_h_tsv.gz` | CAC GWAS (Kavousi et al., Nat Genet 2023) | EBI GWAS Catalog: https://www.ebi.ac.uk/gwas/studies/GCST90278456 |
| `GCST90503074_h_tsv.gz` | SIS GWAS (Gummesson et al., Nat Commun 2025) | EBI GWAS Catalog: https://www.ebi.ac.uk/gwas/studies/GCST90503074 |

### 4. Edit config

Open `R/config.R` and update file paths if needed (defaults assume `data/raw/`).

### 5. Run the pipeline

```r
source("run_pipeline.R")
```

Or run individual steps:

```r
source("R/01_prepare_data.R")
source("R/02_run_mr.R")
# ... etc
```

---

## Methods summary

### Instrument selection
- MiBioGen: p < 1×10⁻⁴, F > 10; 708 SNPs across 66 taxa (MiBioGen-only)
- DMP: p < 1×10⁻⁴, F > 10; 554 SNPs across 65 taxa (DMP-only)
- 37 overlapping genera: pooled after 500 kb distance-based clumping
- Total: 2,020 instruments across 168 taxa

### Outcome meta-analysis
- Fixed-effects IVW meta-analysis of GCST90278456 (CAC score) and GCST90503074 (SIS)
- Allele harmonisation including palindromic SNP handling

### MR estimators
| Method | Assumption | Implementation |
|--------|------------|----------------|
| IVW (primary) | Some pleiotropy allowed (RE) | Weighted regression through origin |
| MR-Egger | Directional pleiotropy detected | WLS with intercept |
| Weighted Median | ≥50% valid instruments | Bootstrap SE (1,000 iterations) |
| Wald Ratio | Single instrument | Direct ratio |

### Multiple testing
- Bonferroni: p < 3.03×10⁻⁴ (0.05 / 165 taxa tested)
- FDR: Benjamini-Hochberg q < 0.05

### Taxonomy audit
- DMP higher-order taxonomy strings (family/order/class/phylum level) excluded
- 36 excluded + 1 unassigned = 37 removed; 81 genus-level taxa retained

### Sensitivity analyses
- Cochran's Q heterogeneity
- MR-Egger intercept (directional pleiotropy)
- Steiger filtering (causal direction)
- Pairwise MVMR within Lachnospiraceae
- Leave-one-out (top 10 genera)
- Funnel plot asymmetry

### Taxonomic pattern analysis
- Family-level mean OR with one-sample t-test
- Group comparison: Lachnospiraceae+Ruminococcaceae vs Other Firmicutes vs Non-Firmicutes
- Mann-Whitney U test (one-sided)
- Chi-square on protective proportions

---

## Output files

| File | Description |
|------|-------------|
| `results/mbg_instruments.csv` | MiBioGen clumped instruments |
| `results/dmp_instruments.csv` | DMP clumped instruments |
| `results/meta_cac_outcome.csv` | Meta-analysed CAC GWAS at instrument positions |
| `results/mr_all_taxa.csv` | MR results for all taxa (pre-audit) |
| `results/mr_clean_genera.csv` | 81 genus-level significant taxa |
| `results/mr_excluded_taxa.csv` | 37 excluded higher-order taxa |
| `results/mr_clean_with_taxonomy.csv` | MR results with phylum/family annotation |
| `results/family_summary.csv` | Family-level mean OR + t-test results |
| `results/taxonomy_stats.rds` | Group comparison test statistics |
| `results/mpgs_weights_plink.txt` | mPGS weights (PLINK-compatible) |
| `results/paf_results.csv` | Population-attributable fractions |
| `results/mpgs_results.rds` | R², AUC, decile OR |
| `figures/fig1_forest_risk.png` | Forest plot: 41 risk-increasing genera |
| `figures/fig2_forest_protective.png` | Forest plot: 40 protective genera |
| `figures/fig_taxonomy.png` | Taxonomic pattern analysis (3-panel) |
| `figures/fig_mpgs_prediction.png` | mPGS prediction (4-panel) |
| `figures/fig_sensitivity.png` | Sensitivity analyses |

---

## Session info

```r
sessionInfo()
# R version 4.3.x
# Key packages: data.table, dplyr, TwoSampleMR, ggplot2, patchwork
```

---

## Citation

If you use this pipeline, please cite:

> [Author Names]. Gut Microbiota Composition and the Risk of Coronary Artery Calcification:
> A Genome-Wide Mendelian Randomization Meta-Analysis. [Journal] [Year].

And the underlying GWAS:
- Kurilshikov A, et al. *Nat Genet.* 2021;53:156–165.
- Lopera-Maya EA, et al. *Nature Genet.* 2022;54:143–151.
- Kavousi M, et al. *Nat Genet.* 2023;55:1651–1664.
- Gummesson A, et al. *Nat Commun.* 2025;16:2266.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contact

Open a GitHub issue or contact [your email] for questions about this pipeline.
