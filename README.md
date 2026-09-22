# Metabolic responses of *Pleurotus djamor* to lignocellulosic substrate gradients

This repository contains the data and R script associated with the study:

**Metabolic responses of *Pleurotus djamor* to lignocellulosic substrate gradients**

## Experimental design

The study evaluated progressive replacement of barley straw with two agro-industrial residues:

- Coconut fiber
- Coffee husk

Replacement levels were 20%, 40%, 60%, 80%, and 100%, with 100% barley straw used as the control.

## Repository contents

### Data/Raw

- `DataBaseGCMS.xlsx`: GC–MS metabolomic dataset before metabolomic filtering and statistical preprocessing.
- `ProductionPdjamor.xlsx`: Productivity data, including yield and biological efficiency.

### Data/Processed

- `MetabolicMatrix.xlsx`: processed metabolomic matrix used for statistical analyses.

### Script

- `MetabolomicAnalysis.R`: R script used for metabolomic data processing and statistical analyses.

## Metabolomic analysis

The workflow includes:

- data quality filtering
- contaminant removal
- aggregation of repeated metabolite annotations
- prevalence filtering
- variance filtering
- Z-score normalization
- PERMANOVA
- betadisper
- PCA
- heatmap analysis
- individual metabolite Kruskal–Wallis tests
- Benjamini–Hochberg FDR correction
- fold-change calculations
- metabolic pathway enrichment analysis

The final metabolomic dataset consisted of 66 metabolites across 30 biological samples.

## Reproducibility

The files provided in this repository are intended to support reproducibility of the analyses reported in the manuscript.
