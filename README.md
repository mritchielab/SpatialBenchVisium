This repository contains code used to perform analysis and figures featured in our paper:

[**Spotlight on 10x Visium: a multi-sample protocol comparison of spatial technologies**](https://www.biorxiv.org/content/10.1101/2024.03.13.584910v1) 

Mei R. M. Du, Changqing Wang, Charity W. Law, Daniela Amann-Zalcenstein, Casey J. A. Anttila, Ling Ling,
Peter F. Hickey, Callum J. Sargeant, Yunshun Chen, Lisa J. Ioannidis, Pradeep Rajasekhar, Raymond Yip, Kelly
L. Rogers, Diana S. Hansen, Rory Bowden, and Matthew E. Ritchie

![Visium data generation and analysis workflow](https://github.com/mritchielab/SpatialBench/blob/main/Visium%20workflow.png) 
Figure created with [BioRender](https://biorender.com).



## Data Availability
Our Visium and 10x scRNA-seq datasets are available from GEO under accession number GSE254652.

Please cite our paper if you use our data and/or scripts in your studies.

## Index

Download the `site` folder and open `index.html` to view HTML reports as a website. Code to produce the reports are stored as Rmarkdown documents in `analysis`. Plots and objects are saved in `output`.

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

Before running the targets project, you will need to replace the `targets_project/data/SpaceRanger` symlink with the appropriate path to the SpaceRanger output folder, and update the `targets_project/data/sample_ids.tsv` file's folder field. 

### Figures

Figures 1-3, Supplementary figures S2-S7: `analysis/figures.R`

Figure 4 & 5, Supplementary figures S8-S11: run `targets::tar_make()` in the `targets_project` folder, figures are saved to `targets_project/output`
