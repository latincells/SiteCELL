# SiteCELL: a protocol to isolate PBMCs directly on remote locations for scRNAseq experiments.  
This repository contains the analysis pipeline for comparing PBMC isolation protocols (SiteCELL vs. Ficoll density gradient) using single-cell RNA-seq data. 
The analysis focuses on cell composition, cell quality metrics, and transcriptomic stress signatures, particularly in samples from diverse ancestry backgrounds. This protocol
is part of the LatinCells project which aims to map an immune atlas of indigenous and admixed Latin American populations.

## How to run

Run scripts in the following order:

### 1. Preprocessing and QC (per dataset)

Each dataset is processed independently, including quality control, filtering, and cell type annotation.

- SiteCELL datasets:
  - `01_preprocessing-SiteCell-Dataset1.R`
  - `02_preprocessing-SiteCell-Dataset2.R`
  - `03_preprocessing-SiteCell-Dataset3.R`

- Public FDG datasets:
  - `04_preprocessing-FDG1-GSE226896.R`
  - `05_preprocessing-FDG2-GSE149689.R`
  - `06_preprocessing-FDG3-GSE161918.R`
  - `07_preprocessing-FDG4-GSE213516.R`
  - `08_preprocessing-FDG5-GSE168732.R`

### 2. Dataset integration

- Merged and integrated datasets 
  - `09_01-all-datasets-merge-and-integration.R`

### 3. Downstream analyses

- Cell quality and composition
  - `09_02-experiment-mitoRNAanalysis.R`
  - `09_03-experiment-cellproportion-calculation.R`

- Stress and functional analyses
  - `09_04-experiment-StressGenesAnalyses.R`
  - `09_05-experiment-subT-diversity.R`

- Replicates and validation
  - `10-experiment-FigS4-replicate-MX0049.R`

- Cross-dataset comparison
  - `11-experiment-comparison-SiteCELL-with-HIHA_dataset.R`

- Sample matching
  - `12_01-preprocessing-matching-samples-with-SiteCELL.R`
  - `12_02-preprocessing-matching-samples-with-FDG.R`

- Supplementary figures
  - `supplementary-fig2.R`



## Methods summary

In this study, we systematically compared two PBMC isolation protocols—LatinCells and Ficoll—using single-cell RNA sequencing (scRNA-seq) data derived from healthy donors with 
diverse ancestry backgrounds. After preprocessing and quality control, individual datasets were processed using a standardized workflow. Cell types were annotated using a consistent reference-based 
approach (Azimuth) to ensure comparability across samples and protocols. To quantify cellular composition, we calculated cell type proportions 
at the individual level of cells recovered per sample. This per-individual approach enabled direct comparison of cellular distributions while accounting for variability in cell quality (genes per cell, UMIs per cell, 
mitochondrial RNA, stress genes) and cell recovery and viability, enabling us to quantify potential artifacts introduced during isolation and processing. This framework allowed 
us to assess whether expression profile biases are protocol-dependent, providing a comprehensive evaluation of protocol performance.


## Data availability

ENA Accession number: E-MTAB-15719 []
Additional files: Zenodo [10.5281/zenodo.19224509]


## Citation

If you use this code, please cite:

Espinosa-Jaime, Aarón, et al. "SiteCELL enables on-site PBMCs purification and cryopreservation for immune single cell profiling of diverse ancestries." bioRxiv (2025): 2025-09. 
