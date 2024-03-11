plot_samples <- function(spes, spot_metric, labs, title,
                         discrete = FALSE, size = 0.25) {
  if (is.list(spes)) {
    p <- sapply(spes,
      function(spe) {
        tibble::tibble(
          counts = spot_metric(spe),
          "x_coor" = spatialCoords(spe)[, 1],
          "y_coor" = spatialCoords(spe)[, 2],
          sample = metadata(spe)$sample
        )
      },
      simplify = FALSE
    ) %>%
      do.call(rbind, .) |>
      ggplot(aes(x = x_coor, y = y_coor, col = counts)) +
      facet_wrap(~sample, ncol = 2)
  } else {
    p <- spatialCoords(spes) |>
      data.frame() |>
      ggplot(aes(
        x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        col = spot_metric(spes)
      )) +
      facet_wrap(~ spes$sample_id, ncol = 2)
  }
  p <- p +
    geom_point(size = size) + coord_fixed() +
    xlab("") + ylab("") + theme_minimal() +
    labs + title +
    if (discrete) {
      scale_colour_discrete()
    } else {
      scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
      )
    }
  return(p)
}

plot_pca_density <- function(spe, dimname = "PCA", dims = c(1, 2), col,
                             density_alpha = 0.3, points_alpha = 0.7) {
  data <- data.frame(
    x = reducedDim(spe, dimname)[, dims[1]],
    y = reducedDim(spe, dimname)[, dims[2]],
    col = if (missing(col)) {
      spe$sample_id
    } else {
      col
    }
  )
  p <- ggplot(data, aes(x = x, y = y, col = col)) +
    geom_point(alpha = points_alpha)
  xdens <- axis_canvas(p, axis = "x") +
    geom_density(data = data, aes(x = x, fill = col, col = col), alpha = density_alpha)
  ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
    geom_density(data = data, aes(x = y, fill = col, col = col), alpha = density_alpha) +
    coord_flip()

  return(
    p |>
      cowplot::insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") |>
      cowplot::insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") |>
      ggdraw()
  )
}

plot_complexheatmap <- function(spe, cluster_label, marker_genes_df,
                                counts_fn = SingleCellExperiment::logcounts) {
  m <- spe[marker_genes_df$EnsembleID, order(cluster_label)] |>
    counts_fn() |>
    as.matrix() |>
    scale()

  col_fun <- circlize::colorRamp2(
    seq(min(m), mean(m) + 2 * sd(m), length.out = 100),
    grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)
  )

  return(ComplexHeatmap::Heatmap(
    m,
    col = col_fun,
    column_split = sort(cluster_label),
    row_split = marker_genes_df$zone,
    cluster_rows = F, cluster_columns = F,
    show_column_names = F,
    use_raster = T, raster_quality = 5,
    raster_resize_mat = mean,
    row_labels = rowData(spe[marker_genes_df$EnsembleID, ])$symbol,
    heatmap_legend_param = list(
      title = "relative expression\n",
      title_gp = grid::gpar(fontsize = 16, fontface = "plain")
    ),
    row_title_rot = 0, column_title_rot = 45
  ))
}

logfc_mtx <- function(spe, cluster_label, marker_genes_df) {
  args.grid <- expand.grid(
    zone = as.character(unique(marker_genes_df$zone)),
    cluster = as.character(unique(cluster_label))
  ) |>
    tibble::as_tibble()

  parallel::mcmapply(function(zone, cluster) {
    rows <- marker_genes_df[marker_genes_df$zone == zone, "EnsembleID"]
    cols <- cluster_label == cluster
    cluster_score <- sum(SingleCellExperiment::counts(spe[rows, cols])) / sum(cols)
    other_score <- sum(SingleCellExperiment::counts(spe[rows, !cols])) / sum(!cols)
    tibble::tibble_row(
      zone = zone,
      cluster = cluster,
      logfc = log2(cluster_score / other_score)
    )
  }, args.grid$zone, args.grid$cluster, mc.cores = 4, SIMPLIFY = F) |>
    dplyr::bind_rows()
}

cowplot_plotVisium <- function(spe, col) {
  spe$col <- col
  cowplot::plot_grid(
    plotlist = unique(spe$sample_id) |>
      lapply(function(i) {
        spe[, spe$sample_id == i] |>
          ggspavis::plotVisium(fill = "col")
      })
  )
}

# https://github.com/XiaoZhangryy/iSC.MEB/blob/233471ca6a96e2b5ada79f150d3574c1cd906e15/R/Visualization.r
iSC.MEB_palette <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
