# file loading
delay_filter <- function(...) {
  args = rlang::enquos(...)
  function(.data) rlang::eval_tidy(dplyr::filter(.data, !!!args))
}

load_samples_tb <- function(sample_ids_file = "data/sample_ids.tsv",
                            filter_func = function(x) {
                              x
                            },
                            sample_name_func = function(x) {
                              gsub("^.*_", "", x$folder)
                            }) {
  folder <- "data/SpaceRanger"
  samples_tb <- readr::read_delim(sample_ids_file) |>
    filter_func()

  samples_tb$raw_spe <-
    sapply(samples_tb$folder, function(x) {
      cat("Loading SCE from", x, "...\n")
      list.files(folder,
        pattern = paste0(x, "$"),
        include.dirs = TRUE,
        recursive = TRUE,
        full.names = TRUE
      ) |>
        SpatialExperiment::read10xVisium(
          type = "sparse", data = "raw",
          images = "hires", load = TRUE
        )
    })
  samples_tb$sample <- sample_name_func(samples_tb)
  names(samples_tb$raw_spe) <- samples_tb$sample
  for (sample in names(samples_tb$raw_spe)) {
    S4Vectors::metadata(samples_tb$raw_spe[[sample]])$sample <- sample
  }
  return(samples_tb)
}

# QC
scuttleFilter <- function(samples_tb,
                          feature_discard_func = function(spe) {
                            rowData(spe)$detected == 0
                          }) {
  samples_tb$raw_spe <- sapply(samples_tb$raw_spe, function(spe) {
    SummarizedExperiment::rowData(spe)$EnsembleID <- rownames(spe)
    # rownames(spe) <- SummarizedExperiment::rowData(spe)$symbol

    # in_tissue
    spe <- spe[, spe$in_tissue]

    # add per-spot QC
    is_mito <- grepl("(^MT-)|(^Mt-)|(^mt-)", SummarizedExperiment::rowData(spe)$symbol)
    spe <- scuttle::addPerCellQC(spe, subsets = list(mito = is_mito))
    qc <- scuttle::perCellQCFilters(spe, sub.fields = "subsets_mito_percent")
    SummarizedExperiment::colData(spe)$discard <- qc$discard

    # add per-feature QC based on kept spots
    spe <- scuttle::addPerFeatureQC(spe[, !colData(spe)$discard])
    rowData(spe)$discard <- feature_discard_func(spe)

    return(spe)
  })

  samples_tb$spe <- sapply(samples_tb$raw_spe, function(spe) {
    spe[!rowData(spe)$discard, !colData(spe)$discard]
  })
  return(samples_tb)
}

runPCA_excludes <- function(spe, exclude_pattern = c(symbol = "^(Ig)|(Gm)"), EnsembleIDs, ...) {
  subset_row <- rep(TRUE, nrow(spe))
  if (!missing(EnsembleIDs)) {
    subset_row <- subset_row & (!rowData(spe)$EnsembleID %in% EnsembleIDs)
  }
  if (!missing(exclude_pattern) && all(!is.na(exclude_pattern))) {
    for (i in names(exclude_pattern)) {
      subset_row <- subset_row & !grepl(exclude_pattern[i], rowData(spe)[, i])
    }
  }
  return(runPCA(spe, subset_row = subset_row, ...))
}

load_filter_seu <- function(samples_tb) {
  stop("WIP")
  samples_tb$seu <-
    sapply(samples_tb$folder, function(x) {
      cat("Loading SEU from", x, "...\n")
      list.files(folder,
        pattern = paste0(x, "$"),
        include.dirs = TRUE,
        recursive = TRUE,
        full.names = TRUE
      ) |>
        SpatialExperiment::read10xVisium(
          type = "sparse", data = "raw",
          images = "lowres", load = TRUE
        )
    })
}

use_mart_sleep <- function(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", max_retry = 5) {
  mart <- NULL
  pause <- 3
  while (is.null(mart) & max_retry > 0) {
    tryCatch(
      expr = {
        mart <- biomaRt::useMart(biomart, dataset = dataset)
      },
      error = function(e) {
        cat(str(e))
        cat("Retry in", pause, "seconds...\n")
        Sys.sleep(pause)
      }
    )
    pause <- pause * 2
    max_retry <- max_retry - 1
  }
  return(mart)
}

add_mart_columns <- function(spe, mart_columns = "entrezgene_id") {
  mart_table <-
    use_mart_sleep() |>
    biomaRt::select(
      keys = rowData(spe)$EnsembleID,
      keytype = "ensembl_gene_id",
      columns = c("ensembl_gene_id", mart_columns)
    )

  rowData(spe) <- rowData(spe) |>
    cbind(mart_table[match(rowData(spe)$EnsembleID, mart_table$ensembl_gene_id), mart_columns, drop = FALSE])
  return(spe)
}

human_to_mouse_genes <- function(df) {
  human <- use_mart_sleep(dataset = "hsapiens_gene_ensembl")
  mouse <- use_mart_sleep(dataset = "mmusculus_gene_ensembl")

  hgnc_symbol <- biomaRt::select(
    human,
    keys = df$human_ensembl_gene_id,
    keytype = "ensembl_gene_id",
    columns = c("ensembl_gene_id", "hgnc_symbol")
  ) |>
    setNames(c("human_ensembl_gene_id", "hgnc_symbol"))

  if (!"ENSG00000073009" %in% hgnc_symbol$human_ensembl_gene_id) {
    hgnc_symbol[nrow(hgnc_symbol) + 1, ] <- c("ENSG00000073009", "IKBKG")
  }

  mgi_symbol <- biomaRt::select(
    mouse,
    keys = hgnc_symbol$hgnc_symbol,
    keytype = "mgi_symbol",
    columns = c("ensembl_gene_id", "mgi_symbol")
  ) |>
    setNames(c("mouse_ensembl_gene_id", "mgi_symbol")) |>
    dplyr::mutate(hgnc_symbol = toupper(mgi_symbol))

  df <- df |>
    merge(hgnc_symbol, all.x = T) |>
    merge(mgi_symbol, all.x = T)
  return(df)
}

to_seu <- function(spe) {
  stopifnot(!anyDuplicated(rownames(spe)))

  s <- SeuratObject::CreateSeuratObject(
    counts = SingleCellExperiment::counts(spe),
    assay = "RNA",
    meta.data = as.data.frame(SingleCellExperiment::colData(spe))
  )
  s@meta.data <- s@meta.data[, c("nCount_RNA", "nFeature_RNA", "in_tissue", "array_row", "array_col")]
  colnames(s@meta.data)[4:5] <- c("row", "col")

  return(s)
}

run_nnSVG <- function(spe_list) {
  sample_ids <- names(spe_list)
  stopifnot(all(sapply(spe_list, \(x) "symbol" %in% colnames(SummarizedExperiment::rowData(x)))))

  res_list <- parallel::mclapply(spe_list, function(spe) {
    SummarizedExperiment::rowData(spe)$gene_name <- SummarizedExperiment::rowData(spe)$symbol
    spe <- nnSVG::filter_genes(
      spe,
      filter_genes_ncounts = 3,
      filter_genes_pcspots = 0.5,
      filter_mito = TRUE
    )
    spe <- spe[, colSums(SingleCellExperiment::counts(spe)) > 0] |>
      scuttle::computeLibraryFactors() |>
      scuttle::logNormCounts() |>
      nnSVG::nnSVG()

    return(rowData(spe))
  }, mc.cores = length(spe_list))

  res_list
}

filter_nnSVG <- function(nnSVG_res_list, padj = 0.05) {
  nnSVG_res_list |>
    lapply(\(x) rownames(x)[x$padj < padj]) |>
    unlist() |>
    unique()
}
