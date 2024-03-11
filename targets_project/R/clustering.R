intersect_spe_rows <- function(samples_tb) {
  rows <- Reduce(intersect, lapply(samples_tb$spe, rownames))
  samples_tb$spe <- sapply(samples_tb$spe, function(x) {
    x[rows, ]
  })
  return(samples_tb)
}

combine_spe <- function(samples_tb) {
  sapply(samples_tb$spe, function(spe) {
    spe$sample_id <- metadata(spe)$sample
    rowData(spe) <- rowData(spe)[, c("symbol", "EnsembleID")]
    colnames(spe) <- paste0(gsub("-1$", "-", colnames(spe)), metadata(spe)$sample)
    return(spe)
  }) %>%
    do.call(cbind, .)
}

get_cluster_label <- function(
    combined_spe, k = 10, resolution_parameter = 0.06,
    use.dimred = "PCA", n_iterations = 100) {
  combined_spe |>
    buildSNNGraph(k = k, use.dimred = use.dimred) |>
    igraph::cluster_leiden(resolution_parameter = resolution_parameter, n_iterations = n_iterations) |>
    igraph::membership() |>
    factor()
}

run_iSC_MEB <- function(
    seu_list, customGenelist, k, n_pca = 10, maxIter = 20) {
  iSC.MEB::CreateiSCMEBObject(seu_list, customGenelist = customGenelist) |>
    iSC.MEB::CreateNeighbors(platform = "Visium") |>
    iSC.MEB::runPCA(npcs = n_pca, pca.method = "APCA") |>
    iSC.MEB::SetModelParameters(maxIter = maxIter, coreNum = length(seu_list)) |>
    iSC.MEB::iSCMEB(K = k)
}
