# TNBC-scRNAseq-CNV-Analysis
Computational analysis pipeline for single-cell RNA-sequencing data from human TNBC tumors, including CNV inference, tumor-cell annotation, NK and tumor population analysis, immune-cell analysis, and differential expression.

Reproducible R workflow for single-cell RNA-seq, CNV, tumor-cell refinement, CSC/FZD7 analysis, WNT-associated gene analysis, SOCS3 analysis, and NK-cell analysis in triple-negative breast cancer (TNBC).

## Overview

This repository contains the computational scripts used for the analysis of
single-cell RNA-seq data from TNBC samples.

The workflow includes:

1. Raw 10X Genomics data processing
2. Ambient RNA correction using SoupX
3. Initial and sample-specific quality control
4. Seurat normalization and dimensionality reduction
5. CopyKAT-based CNV inference
6. Cell-type annotation
7. Refinement of the aneuploid tumor-cell population
8. CSC and FZD7 classification
9. WNT-associated gene expression analysis
10. SOCS3 expression analysis
11. NK-cell identification, refinement, and reclustering
12. Response-wise statistical analysis using Wilcoxon rank-sum tests
13. Multiple-testing correction using the Benjamini-Hochberg method

## Repository structure

```text
TNBC-scRNAseq-CNV-Analysis/
├── README.md
├── LICENSE
├── .gitignore
├── metadata/
│   └── sample_metadata.csv
└── scripts/
    └── FINAL_SC_RNA_CNV_ANALYSIS.R
