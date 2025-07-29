# SpatialBenchVisium

This repository contains code used to perform analysis and create figures featured in our paper:

[Du, Wang, Law, Amann-Zalcenstein *et al.* (2025) **Benchmarking spatial transcriptomics technologies with the multi-sample SpatialBenchVisium dataset**, Genome Biol 26:77.](https://doi.org/10.1186/s13059-025-03543-4)

![Visium data generation and analysis workflow](https://github.com/mritchielab/SpatialBench/blob/main/Visium%20workflow.png) 
Figure created with [BioRender](https://biorender.com).

## Data Availability

Our processed Visium and 10x scRNA-seq datasets, along with the code are available from zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12788291.svg)](https://doi.org/10.5281/zenodo.12788291), data is also accessible through GEO: [GSE254652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254652).

Please cite our [paper](https://doi.org/10.1186/s13059-025-03543-4) if you use our data and/or scripts in your research.

## Index

Code to produce the reports are stored as Rmarkdown documents in `analysis`. Objects are saved in `output`.

### Pre-processing

#### Spatial
Sample 709 FFPE CA: `analysis/EDA_709_FFPE_CA.Rmd`

Sample 713 FFPE CA: `analysis/EDA_713_FFPE_CA.Rmd`

#### scRNA-seq
`analysis/sc_preprocessing.Rmd`

### Downstream analysis

Multi-sample feature selection, clustering, cell type deconvolution: `analysis/FFPE_CA_multi-sample.Rmd`

Pseudo-bulk differential expression analysis: as [targets](https://docs.ropensci.org/targets/) project under the `targets_project` folder

#### Running targets

Simply navigate to the `targets_project` folder from the zenodo tarball `SpatialBenchVisium.tar.gz` and run `targets::tar_make()` to run the entire pipeline, outputs are saved in the `targets_project/output` folder.

### Figures

Figures 1-3, Supplementary figures S2-S7: `analysis/figures.R`

Figure 4 & 5, Supplementary figures S8-S11: run `targets::tar_make()` in the `targets_project` folder, figures are saved to `targets_project/output`
