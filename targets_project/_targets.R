library(targets)
library(tarchetypes)
tar_option_set(
  packages = c(
    "tibble", "SpatialExperiment", "RColorBrewer", "gridExtra",
    "grDevices", "SingleCellExperiment", "ggplot2", "tidyverse",
    "scater", "scran", "ggspavis", "limma", "edgeR", "cowplot",
    "patchwork", "tweeDEseqCountData"
  ),
  format = "qs",
  # memory = "transient", garbage_collection = TRUE,
  storage = "worker", retrieval = "worker"
)
tar_source()

unlist(
  lapply(
    X = list.files(
      path = file.path("targets"),
      pattern = "\\.R$",
      full.names = TRUE,
      all.files = TRUE,
      recursive = TRUE
    ),
    FUN = function(path) {
      eval(
        expr = parse(text = readLines(con = path, warn = FALSE)),
        envir = targets::tar_option_get(name = "envir")
      )
    }
  )
)
