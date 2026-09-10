# TNBC-scRNAseq-CNV-Analysis
Computational analysis pipeline for single-cell RNA-sequencing data from human TNBC tumors, including CNV inference, tumor-cell annotation, NK and tumor population analysis, immune-cell analysis, and differential expression.

Reproducible R workflow for single-cell RNA-seq, CNV, tumor-cell refinement, CSC/FZD7 analysis, WNT-associated gene analysis, SOCS3 analysis, and NK-cell analysis in triple-negative breast cancer (TNBC).

## Overview

This repository contains the computational scripts used for the analysis of
single-cell RNA-seq data from TNBC samples.

The workflow includes:

Raw 10X Genomics data processing
Ambient RNA correction using SoupX
Initial and sample-specific quality control
Seurat normalization and dimensionality reduction
CopyKAT-based CNV inference
Cell-type annotation
Refinement of the aneuploid tumor-cell population
CSC and FZD7 classification
WNT-associated gene expression analysis
SOCS3 expression analysis
NK-cell identification, refinement, and reclustering
Response-wise statistical analysis using Wilcoxon rank-sum tests
Multiple-testing correction using the Benjamini-Hochberg method

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
