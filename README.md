# Development of a Circulating Protein–Metabolite Signaling Network

## Overview

This project develops a statistically robust framework to identify functional relationships between circulating proteins and metabolites in human blood using population-level data.

Modern biomedical research increasingly emphasizes **integrative multi-omics analysis**. While proteomics and metabolomics are often studied independently, this project aims to construct a **cross-omic network** that captures how these molecular layers co-vary.

* **Proteins** reflect biological regulation ("intent")
* **Metabolites** reflect biochemical outcomes ("effect")

The resulting network provides a systems-level view of human physiology and serves as a foundation for future disease-focused analyses.

---

## Data Description

* ~10,000 individuals (after quality control)
* ~2,900 proteins (Olink platform)
* ~249 metabolites (NMR platform)

Data preprocessing includes:

* Missingness assessment
* Outlier detection
* Covariate adjustment (age, sex)

---

## Project Pipeline

The analysis is structured into reproducible blocks:

### Block 0 — Project Setup

* R project initialization
* Reproducibility setup using `renv`
* Organized directory structure

### Block 1 — Quality Control (QC)

* Feature-wise and sample-wise missingness analysis
* Visualization of missing data patterns

### Block 2 — Imputation

* Matrix completion using `softImpute`
* Efficient handling of high-dimensional proteomics data

### Block 3 — Outlier Detection

* PCA-based dimensionality reduction
* Mahalanobis distance for multivariate outlier removal

### Block 4 — Covariate Adjustment

* Linear regression per feature
* Adjustment for:

  * Sex
  * Age (natural splines)

### Block 5 — Protein–Metabolite Associations

* Pearson and Spearman correlations
* Multiple testing correction (FDR)
* Thresholding based on effect size and significance

### Block 6 — Multivariate Integration (CCA)

* Ridge-regularized CCA
* Classical CCA (comparison)
* Sparse CCA (interpretability)

### Block 7 — Network Construction

* Bipartite protein–metabolite network
* Community detection (Louvain method)
* Functional enrichment (GO, KEGG)

---

## Repository Structure

```
Thesis_R_Project/
│
├── data/            # Raw and intermediate data (excluded from GitHub)
├── scripts/         # Main analysis scripts (pipeline blocks)
├── R/               # Helper functions
├── results/         # Output tables and intermediate results
├── figures/         # Generated figures for analysis and reports
├── reports/         # R Markdown / thesis documents
├── logs/            # Execution logs
├── renv/            # Reproducible environment
└── deliverables/    # Final outputs and exports
```

---

## Reproducibility

This project uses:

* **R**
* **renv** for dependency management

To reproduce the environment:

```r
install.packages("renv")
renv::restore()
```

---

## Key Methods

* Pearson correlation
* False Discovery Rate (Benjamini–Hochberg)
* PCA + Mahalanobis distance
* Canonical Correlation Analysis (CCA)
* Network analysis (Louvain clustering)
* Functional enrichment (GO, KEGG)

---

## Results Summary

* Large-scale protein–metabolite association network constructed
* Identification of biologically meaningful modules
* Strong enrichment of immune and metabolic pathways
* Evidence of coordinated regulation between proteins and metabolites

---

## Notes

* Raw data is not included due to size and privacy constraints
* Only processed outputs and scripts are provided

---

## Author

Sangeeth S Menon
MSc Statistics and Machine Learning
Linköping University

---

## Future Work

* Robustness analysis (Block 8)
* Disease-specific module analysis
* Extension to longitudinal data
* Integration with additional omics layers

---
