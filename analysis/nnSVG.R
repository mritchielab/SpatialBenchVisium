# Multi-sample feature selection with nnSVG
# Author: Mei Du

library(here)
library(SpatialExperiment)
library(Seurat)
library(scuttle)
# remotes::install_github("lmweber/nnSVG")
# parallelisation
library(foreach)
library(doParallel)

# set number of cores
workers <- 4
cl <- parallel::makeCluster(workers)
# register the cluster for using foreach
doParallel::registerDoParallel(cl)

# Load data
# pre-processed SpatialExperiment objects
spe_709 <- readRDS(here("analysis","output","RDS","spe_709_FFPE_CA_norm.RDS"))
spe_713 <- readRDS(here("analysis","output","RDS","spe_713_FFPE_CA_norm.RDS"))

# Combining samples with SpatialExperiment's cbind()
spe_709$sample_id <- "709"
spe_713$sample_id <- "713"
spe <- cbind(spe_709, spe_713)
table(colData(spe)$sample_id)

# HVGs
hvg_709 <- read.csv(here("analysis","output","HVGs","709_FFPE_CA_top_HVGs.csv"))
hvg_713 <- read.csv(here("analysis","output","HVGs","713_FFPE_CA_top_HVGs.csv"))

# nnSVG
library(nnSVG)

# Run nnSVG once per sample and store lists of top SVGs
sample_ids <- unique(colData(spe)$sample_id)
colnames(rowData(spe)) <- c("gene_name","Chr")

res_list <- foreach(s=seq_along(sample_ids), .multicombine=TRUE, .packages=c("SpatialExperiment","Seurat","scuttle","nnSVG")) %dopar% {
  # select sample
  ix <- colData(spe)$sample_id == sample_ids[s]
  spe_sub <- spe[, ix]

  # run nnSVG filtering for mitochondrial genes and low-expressed genes
  spe_sub <- filter_genes(
    spe_sub,
    filter_genes_ncounts = 3,
    filter_genes_pcspots = 0.5,
    filter_mito = TRUE
  )

  # remove any zeros introduced by filtering
  ix_zeros <- colSums(counts(spe_sub)) == 0
  if (sum(ix_zeros) > 0) {
    spe_sub <- spe_sub[, !ix_zeros]
  }

  dim(spe_sub)

  # re-calculate logcounts after filtering
  spe_sub <- computeLibraryFactors(spe_sub)
  spe_sub <- logNormCounts(spe_sub)


  # run nnSVG
  set.seed(123)
  spe_sub <- nnSVG(spe_sub)

  # store results for this sample
  rowData(spe_sub)
}

saveRDS(res_list, here("output","rds","nnSVG_FFPE_CA.RDS"))
