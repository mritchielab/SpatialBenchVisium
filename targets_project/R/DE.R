aggregate_sample_bulk_counts <- function(combined_spe, set.names = "symbol", symbols_pattern = "^Ig") {
  combined_spe <- combined_spe[!is.na(rowData(combined_spe)[, set.names]), ]
  if (!is.na(symbols_pattern)) {
    combined_spe <- combined_spe[!grepl(symbols_pattern, rowData(combined_spe)$symbol), ]
  }
  rownames(combined_spe) <- rowData(combined_spe)[, set.names]
  sample_bulk_counts <- combined_spe |>
    counts() |>
    t() |>
    as.matrix() |>
    by(combined_spe$sample_id, colSums) |>
    sapply(identity) |>
    by(rownames(combined_spe), colSums) |>
    sapply(identity) |>
    t()
  return(sample_bulk_counts)
}
