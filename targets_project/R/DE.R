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

my_enrichGO <- function(gene_list, key_type = "ensembl_gene_id", ont = "BP") {
  entriz_ids <- use_mart_sleep(dataset = "mmusculus_gene_ensembl") |>
    biomaRt::select(
      keys = gene_list,
      keytype = key_type,
      columns = c(key_type, "entrezgene_id")
    ) |>
    na.omit()
  clusterProfiler::enrichGO(
    gene = as.character(entriz_ids$entrezgene_id),
    OrgDb = org.Mm.eg.db::org.Mm.eg.db, # mouse
    keyType = "ENTREZID",
    ont = ont,
    readable = TRUE
  )
}

my_enrichKEGG <- function(gene_list, key_type = "ensembl_gene_id") {
  entriz_ids <- use_mart_sleep(dataset = "mmusculus_gene_ensembl") |>
    biomaRt::select(
      keys = gene_list,
      keytype = key_type,
      columns = c(key_type, "entrezgene_id")
    ) |>
    na.omit()
  clusterProfiler::enrichKEGG(
    gene = as.character(entriz_ids$entrezgene_id),
    organism = "mmu",
    pvalueCutoff = 0.05
  )
}
