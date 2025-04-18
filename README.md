## Introduction
This repository contains scripts for running ATAC-seq downstream analysis. It includes scripts for setting up the project environment, performing differential analysis, generating figures, and creating tables. The code is organized to ensure a clear workflow, where figures depend on both analysis and tables.

---

## Overview
To preprocess the data from each assay, we used the ATAC-seq (https://nf-co.re/atacseq/2.1.2/) pipeline. Then we performed downstream analysis using the processed data. For ATAC-seq, we used Diffbind (v 3.12.0) (using DESeq2 for differential analysis), ChIPSeeker (v 1.38.0) for differential peak annotation, pheatmap (v 1.012) for generating heatmaps and clusterProfiler (v 4.10.10) for functional enrichment (Gene Ontology, Reactome PA, KEGG). 

---

## Table of Contents

- [Setup](#setup)
- [Preprocessing](#preprocessing)
- [Data](#data)
- [Analysis](#analysis)
- [Figures](#figures)
- [Tables](#tables)
- [RObjects](#robjects)


---

## Setup
The majority of this project was done using R (v 4.3.2) and some shell scripting. 

Before running any analysis or generating figures, you need to set up the environment using the `scripts/setup/set_up_environment.R` script. This will ensure all dependencies are installed. This will also create all required directories.

### To run the setup:

1. Clone the repository:
   ```bash
   git clone https://github.com/luulintran/ATAC_seq.git
   cd ATAC_seq
   ```

2. Run the setup script to install the required dependencies:
    ```
    Rscript "scripts/setup/set_up_environment.R"
    ```

## Preprocessing

The `scripts/preprocess/` directory contains scripts to run the Nf-core pipelines for preprocessing ATAC-seq data, as well as filtering fragments for ATAC-seq.

Note: You will need to include a config file to run the Nf-core pipeline scripts. 

## Data

Processed data (bam files and peak files) should be put in the `data/processed_data` directory. The `SampleSheet.csv` should be put in the `data/meta_data` directory. 

## Analysis

The `scripts/analysis/` directory contains scripts for performing downstream analysis using Diffbind, ChIPseeker, clusterProfiler, HOMER, etc.

Within each `analysis` subdirectory, for example `analysis/diffbind`, there will be `config.R` scripts where you can change file paths for input files and output directories, as well as variables, such as control_group and mutant_group.

`scripts/analysis/` contains scripts for differential peak analysis (`diffbind`) and functional enrichment analysis such as Gene Ontology Enrichment Analysis, Reactome Pathway Analysis, and KEGG Pathway Analysis (`functional_enrichment`). These are all done using R. The final directory `homer` contains scripts for Homer motif enrichment analysis using shell scripting.

## Figures

The `results/figures/` directory will contain .png figures that were generated by the scripts in `scripts/analysis/`.

## Tables

The `results/tables/` directory will contain tables (bed files, csv files, txt files) that were generated by the scripts in `scripts/analysis/`.

## RObjects

The `results/r_objects/` directory will contain R objects storing data from Diffbind analysis, which can be used to generate figures and tables.

## Notes
Ensure you have the appropriate data files in `data/processed_data/`. And look at the DiffBind vignette to make sure your `SampleSheet.csv` file is properly formatted

