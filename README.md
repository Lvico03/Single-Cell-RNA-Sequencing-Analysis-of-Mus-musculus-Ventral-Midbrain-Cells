# Single-Cell RNA-Seq Analysis of *Mus musculus* Ventral Midbrain

Single-cell RNA sequencing (scRNA-seq) analysis of adolescent mouse ventral
midbrain, carried out for the Genomics & Transcriptomics course of the MSc in
Bioinformatics for Computational Genomics (Università degli Studi di Milano /
Politecnico di Milano).

## Overview

This project analyses a 10x Chromium scRNA-seq dataset from the ventral midbrain
of an adolescent (P20–30) *Mus musculus*, following the standard Seurat workflow:
quality control, normalization, dimensionality reduction, clustering, marker-gene
identification and manual cell-type annotation. The goal is to recover the
biological cell types present in the tissue and compare them against reference
annotations.

## Dataset

- **Source:** [PanglaoDB](https://panglaodb.se/)
- **SRA:** SRA670243 &nbsp;|&nbsp; **SRS:** SRS3078084 &nbsp;|&nbsp; **SRR:** SRR6877265
- **Tissue:** Ventral midbrain
- **Organism:** *Mus musculus* (adolescent)
- **Protocol:** 10x Chromium &nbsp;|&nbsp; **Instrument:** Illumina HiSeq 2500
- **Original study:** Zeisel et al. (2018), *Molecular Architecture of the Mouse
  Nervous System*, *Cell* 174(4):999–1014.
  [DOI: 10.1016/j.cell.2018.06.021](https://doi.org/10.1016/j.cell.2018.06.021)

## Workflow

1. **Loading & preprocessing** — count matrix loaded from the PanglaoDB
   `.RData` object; gene names cleaned (Ensembl IDs stripped, duplicates made
   unique); Seurat object created (`min.cells = 3`, `min.features = 200`).
2. **Quality control** — filtering on number of detected genes
   (200 < nFeature_RNA < 2500) and mitochondrial content (percent.mt < 5%).
   Mitochondrial and ribosomal fractions computed and inspected.
3. **Normalization** — log-normalization (scale factor 10,000).
4. **Feature selection & scaling** — 2,000 most variable genes; data scaled.
5. **Dimensionality reduction** — PCA; number of PCs chosen using the 75%
   cumulative-variance rule (10 PCs) and the elbow plot (15 PCs).
6. **Clustering** — shared-nearest-neighbour graph + Louvain clustering, tested
   at 10 and 15 PCs (resolution 0.5); UMAP and t-SNE for visualization.
7. **Bias assessment** — clusters checked for mitochondrial, ribosomal and
   cell-cycle artefacts.
8. **Marker genes & annotation** — `FindAllMarkers` (one-vs-all, Wilcoxon test),
   pairwise comparisons for related clusters, manual cell-type assignment.

## Key results

- **6,478 cells** retained after QC (from 9,595 raw / 7,965 after object creation).
- **18 clusters** (15-PC solution) annotated into the major ventral-midbrain
  cell types.
- Cell types recovered: **astrocytes, microglia, macrophages, endothelial cells,
  pericytes, smooth muscle cells, fibroblasts, ependymal cells, the full
  oligodendrocyte lineage** (OPCs → committed OPCs → newly formed → mature →
  myelinating) **and neurons** (dopaminergic and GABAergic).
- A continuous **oligodendrocyte differentiation trajectory** is visible as an
  arc in the UMAP, confirmed by stage-specific markers.
- **Dopaminergic neurons** (Th+, Snca+) — the signature ventral-midbrain
  population affected in Parkinson's disease — were identified.
- No cluster was driven by technical bias (mitochondrial / ribosomal /
  cell-cycle).
  <img width="3000" height="2400" alt="image" src="https://github.com/user-attachments/assets/ced4ea22-8f94-405c-a7fb-56ce5396f1a5" />


## Tools

- **R** with **Seurat** (v5)
- `dplyr`, `ggplot2`, `patchwork`, `stringr`

## Repository contents

- `script` — sc_rna_sequencing_script 
- `plots/` — output plots (QC, PCA, clustering, marker expression, final UMAP)
- `README.md` — ReadMe
- `SRA670243_SRS3078084.sparse.RData` — Data

## Author

Ludovico Salmin — MSc Bioinformatics for Computational Genomics
