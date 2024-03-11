plot_splsda <- function(splsda_res, spe, label, rows.x = 1:5, rows.y = 6:10, size = 0.5, alpha = 0.8) {
  weights <- mixOmics::selectVar(splsda_res)
  spe <- spe[match(rownames(weights$value)[c(rows.x, rows.y)], rowData(spe)$EnsembleID), ]
  data.frame(
    x = (weights$value[rows.x, ] %*% logcounts(spe[seq_along(rows.x), ]))[1, ],
    y = (weights$value[rows.y, ] %*% logcounts(spe[seq(from = length(rows.x) + 1, along.with = rows.y), ]))[1, ],
    col = label
  ) |>
    ggplot(aes(x = x, y = y, col = col)) +
    geom_point(size = size, alpha = alpha) -> p1

  spe |>
    rowData() |>
    cbind(weights$value[c(rows.x, rows.y), , drop = F]) |>
    as_tibble() |>
    tableGrob() -> p2

  return(cowplot::plot_grid(p1, p2, ncol = 2))
}
