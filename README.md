# RNA-seq Expression QC and Clustering Pipeline

## Overview

This project simulates gene expression data and implements a modular
analysis pipeline including:

- Data simulation
- Missingness injection
- Normalization (CPM + log transform)
- Quality control metrics
- Principal Component Analysis (PCA)
- k-means clustering
- Visualization of clustering structure

The goal is to demonstrate foundational bioinformatics data processing
and unsupervised analysis workflows in R.

---

## Project Structure

scripts/
- simulate_data.R
- preprocessing.R
- qc_analysis.R
- dim_reduction.R
- run_pipeline.R

report/
- analysis_report.Rmd

---

## Methods

### 1. Simulation
Gene expression values are simulated using a log-normal distribution.

### 2. Missingness
Missing values are introduced randomly at a specified fraction.

### 3. Normalization
Counts are scaled to counts-per-million (CPM) and log-transformed.

### 4. QC Metrics
Per-sample and per-gene metrics include:
- Missing fraction
- Mean/median expression
- Detection rates
- Correlation to mean sample

### 5. Dimensionality Reduction
PCA is applied to normalized expression data after mean imputation.

### 6. Clustering
k-means clustering is performed on selected principal components.

---

## How to Run

```r
source("R/run_pipeline.R")

res <- run_pipeline(
  num_genes = 2000,
  num_samples = 30,
  missing_frac = 0.1,
  seed = 1
)
