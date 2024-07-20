list(
  tar_target(
    spes_koctrl_tb,
    all_samples |>
      filter(experiment == "OCT Experiment 2" & group != "WT OCT") |>
      mutate(sample = gsub(".*_", "", folder)) |>
      scuttleFilter(feature_discard_func = function(x) {
        FALSE
      })
  ),
  tar_target(
    spe_koctrl,
    spes_koctrl_tb |>
      intersect_spe_rows() |> # subset to intersect of genes
      combine_spe() |>
      # feature QC
      (function(spe) {
        spe$sample_id <- gsub("(A|B)$", "", spe$sample_id)
        spe <-
          scuttle::addPerCellQC(spe, subsets = list(
            mito = grepl("(^MT-)|(^Mt-)|(^mt-)", SummarizedExperiment::rowData(spe)$symbol),
            rps = grepl("^R|rps|l", SummarizedExperiment::rowData(spe)$symbol)
          )) |>
          scuttle::addPerFeatureQC(
            subsets =
              sapply(levels(factor(spe$sample_id)), function(x) {
                x == factor(spe$sample_id)
              }, simplify = F)
          )
        return(spe)
      })() |>
      # raw spe
      computeLibraryFactors() |>
      logNormCounts() |>
      add_mart_columns(mart_columns = "chromosome_name") |>
      (\(spe) {
        genes <-
          rowData(spe) |>
          as_tibble() |>
          # at least 2% spots detection in at least 4 samples
          mutate(low_expr = rowSums(across(matches("subsets_.*_detected")) > 2) < 4) |>
          mutate(excluded = chromosome_name %in% c("X", "Y") | grepl("^Ig", symbol)) |>
          mutate(genes = !(low_expr | excluded)) |>
          pull(genes)
        return(runPCA(spe, subset_row = genes))
      }())() |>
      runUMAP(dimred = "PCA") |>
      (\(x) { # add group labels
        x$group <- factor(ifelse(x$sample_id %in% c("167", "168"), "ko", "ctrl"))
        return(x)
      })()
  ),
  tar_target(iSCMEB_koctrl, 
    {
      set.seed(123)
      if (file.exists('output/iSCMEB_koctrl.qs')) {
        qs::qread('output/iSCMEB_koctrl.qs')
      } else {
        spes_koctrl_tb$spe |>
          sapply(to_seu) |>
          run_iSC_MEB(customGenelist = filter_nnSVG(run_nnSVG(spes_koctrl_tb$spe)), k = 8)
      }
    }
  ),
  tar_target(
    cluster_koctrl,
    iSCMEB_koctrl@resList@idents |>
      unlist() |>
      factor()
  ),
  tar_target(
    plot_umap_koctrl,
    format = "file",
    {
      p <- iSCMEB_koctrl |>
        iSC.MEB::CalculateUMAP(reduction = "iSCMEB", n_comp = 2) |>
        iSC.MEB::LowEmbedPlot(item = "cluster", reduction = "UMAP2", point_size = 0.1, point_alpha = 0.8) +
        theme(
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20)
        ) +
        guides(color = guide_legend(override.aes = list(size = 4)))
      ggsave(file.path("output", "KO_CTL_UMAP.pdf"), p, width = 12, height = 12, dpi = 600, units = "cm")
      file.path("output", "KO_CTL_UMAP.pdf")
    }
  ),
  tar_target(
    plot_cluster_koctrl,
    format = "file",
    {
      tb <- spe_koctrl |>
        spatialCoords() |>
        data.frame() |>
        as_tibble() |>
        mutate(
          group = spe_koctrl$group,
          sample_id = spe_koctrl$sample_id,
          cluster = cluster_koctrl
        )

      theme <- theme_minimal() +
        theme(
          plot.margin = unit(c(1, 0, 0, 0), "cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
        )

      p1 <- tb |>
        filter(group == "ko") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.1) +
        coord_fixed() +
        facet_wrap(~sample_id, ncol = 1) +
        theme +
        scale_color_manual(values = iSC.MEB_palette(length(levels(cluster_koctrl)))) +
        xlab("") +
        ylab("") +
        ggtitle("KO")
      p2 <- tb |>
        filter(group == "ctrl") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.1) +
        coord_fixed() +
        facet_wrap(~sample_id, ncol = 1) +
        theme +
        scale_color_manual(values = iSC.MEB_palette(length(levels(cluster_koctrl)))) +
        xlab("") +
        ylab("") +
        guides(color = guide_legend(override.aes = list(size = 4))) +
        ggtitle("CTL")

      p <- (p1 + guides(color = "none") | p2) + plot_layout(guides = "collect")
      ggsave(file.path("output", "KO_CTL_cluster_plot.pdf"), p, width = 15, height = 10, dpi = 600, units = "cm")
      file.path("output", "KO_CTL_cluster_plot.pdf")
    }
  ),
  tar_target(
    plot_zone_koctrl,
    format = "file",
    {
      tb <- spe_koctrl |>
        spatialCoords() |>
        data.frame() |>
        as_tibble() |>
        mutate(
          group = spe_koctrl$group,
          sample_id = spe_koctrl$sample_id,
          cluster = zone_koctrl
        )

      theme <- theme_minimal() +
        theme(
          plot.margin = unit(c(1, 0, 0, 0), "cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
        )

      p1 <- tb |>
        filter(group == "ko") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.1) +
        coord_fixed() +
        facet_wrap(~sample_id, ncol = 1) +
        theme +
        scale_color_manual(values = mei_palette) +
        guides(color = guide_legend(override.aes = list(size = 4))) +
        xlab("") +
        ylab("") +
        ggtitle("KO")
      p2 <- tb |>
        filter(group == "ctrl") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.1) +
        coord_fixed() +
        facet_wrap(~sample_id, ncol = 1) +
        theme +
        scale_color_manual(values = mei_palette) +
        guides(color = guide_legend(override.aes = list(size = 4))) +
        xlab("") +
        ylab("") +
        ggtitle("CTL")

      p <- (p1 + guides(color = "none") | p2) + plot_layout(guides = "collect")
      ggsave(file.path("output", "KO_CTL_zone_plot.pdf"), p, width = 15, height = 10, dpi = 600, units = "cm")
      file.path("output", "KO_CTL_zone_plot.pdf")
    }
  ),
  tar_target(
    heatmap_koctrl,
    format = "file",
    {
      p <- logfc_mtx(spe_koctrl, cluster_koctrl, zone_markers_df) |>
        ggplot(aes(
          x = factor(cluster, levels = levels(cluster_koctrl)),
          y = factor(zone, levels = rev(levels(zone_markers_df$zone))), fill = logfc
        )) +
        geom_tile(
          color = "white",
          lwd = 0.5,
          linetype = 1
        ) +
        theme_void() +
        theme(
          axis.title = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12), legend.title = element_blank(),
          legend.text = element_text(size = 12)
        ) +
        geom_text(aes(label = format(round(logfc, digits = 1), nsmall = 1)), color = "black", size = 4) +
        scale_fill_distiller(palette = "RdBu") +
        coord_fixed()
      ggsave(file.path("output", "KO_CTL_heatmap.pdf"), p, width = 15, height = 15, dpi = 600, units = "cm")
      file.path("output", "KO_CTL_heatmap.pdf")
    }
  ),
  # manual annotation for cluster
  tar_target(zone_koctrl, cluster_koctrl |>
    (\(x) {
      x$zone <- factor(
        case_when(
          x == 1 ~ "B cell",
          x == 2 ~ "Erythrocyte",
          x == 3 ~ "Neutrophil",
          x == 4 ~ "Erythrocyte",
          x == 5 ~ "Plasma cell",
          x == 6 ~ "Erythrocyte",
          x == 7 ~ "T cell",
          x == 8 ~ "Germinal centre"
        )
      )
    })()),
  ##
  # pseudo-bulk DE
  tar_target(
    puseudo_bulk_counts_koctrl,
    sapply(levels(zone_koctrl),
      function(x) {
        aggregate_sample_bulk_counts(spe_koctrl[, zone_koctrl == x], set.names = "EnsembleID")
      },
      simplify = F
    )
  ),
  tar_target(
    vooms_koctrl,
    sapply(puseudo_bulk_counts_koctrl,
      function(x) {
        dgelist <- DGEList(
          counts = x,
          samples = colnames(x),
          gene = rowData(spe_koctrl)[
            match(rownames(x), rowData(spe_koctrl)$EnsembleID),
            c("symbol", "chromosome_name")
          ],
          group = factor(ifelse(colnames(x) %in% c("167", "168"), "ko", "ctrl"))
        )
        # dgelist <- dgelist[filterByExpr(dgelist), , keep.lib.sizes = FALSE]
        dgelist <- dgelist[
          filterByExpr(dgelist), , # bias?
          keep.lib.sizes = FALSE
        ]
        dgelist <- calcNormFactors(dgelist, method = "TMM")
        # means model
        design <- model.matrix(~ 0 + dgelist$samples$group, data = dgelist$samples)
        colnames(design) <- gsub(".*\\$", "", colnames(design)) |>
          gsub(" .*$", "", x = _)
        dge_v <- voom(dgelist, design, save.plot = T, plot = F, span = 0.2)
        return(dge_v)
      },
      simplify = F
    )
  ),
  tar_target(
    efits_koctrl,
    sapply(names(vooms_koctrl),
      function(x) {
        efit <- lmFit(vooms_koctrl[[x]], vooms_koctrl[[x]]$design) |>
          contrasts.fit(contrasts = makeContrasts(
            Kovsctrl = "groupko - groupctrl",
            levels = vooms_koctrl[[x]]$design
          )) |>
          eBayes()
        efit$voom.line <- vooms_koctrl[[x]]$voom.line
        efit$voom.xy <- vooms_koctrl[[x]]$voom.xy
        efit$cluster <- x
        return(efit)
      },
      simplify = F
    )
  ),
  tar_target(
    prev_de,
    read_csv("data/GCB-DE-D15.TbetKO_GCB.v.WT_GCB.csv") |>
      filter(adj.P.Val < 0.05)
  ),
  tar_target(
    roasts_koctrl,
    sapply(levels(zone_koctrl),
      function(x) {
        roast(vooms_koctrl[[x]],
          design = efits_koctrl[[x]]$design,
          index = efits_koctrl[[x]]$gene %>%
            mutate(idx = seq.int(nrow(.))) %>%
            filter(symbol %in% prev_de$Symbol) %>%
            pull(idx),
          gene.weights = efits_koctrl[[x]]$gene %>%
            filter(symbol %in% prev_de$Symbol) %>%
            left_join(prev_de, by = c("symbol" = "Symbol")) %>%
            pull(t),
          contrast = efits_koctrl[[x]]$contrasts,
          nrot = 99999
        )
      },
      simplify = F
    )
  ),
  tar_target(
    roast_table_koctrl,
    sapply(roasts_koctrl, function(x) {
      x$p.value["Up", ]
    }) |>
      t() |>
      cbind(Number.of.genes = sapply(roasts_koctrl, function(x) {
        x$ngenes.in.set
      }))
  ),
  tar_target(koctrl_DE_genes,
    format = "file",
    (\(efit) {
      file_name <- paste0(file.path("output", "koctrl_"), gsub(" ", "_", efit$cluster), "_DE_genes.csv")
      limma::topTable(efit, number = Inf, p.value = 0.05) |>
        mutate(Overlap_Ly_et_al = symbol %in% prev_de$Symbol) |>
        write_csv(file_name)
      return(file_name)
    })(efits_koctrl[[1]]),
    pattern = map(efits_koctrl),
    iteration = "list"
  ),
  tar_target(KO_CTL_MAplot_Bcell,
    format = "file",
    efits_koctrl[["B cell"]] |>
      (\(efit) {
        pdf(file.path("output", "KO_CTL_MAplot_Bcell.pdf"), width = 6, height = 5)
        status <- efit$genes |>
          dplyr::left_join(prev_de, by = c("symbol" = "Symbol")) |>
          mutate(status = case_when(
            is.na(logFC) ~ "Other",
            symbol == "Tbx21" ~ "T-bet",
            logFC > 0 ~ "T-bet up-regulated",
            logFC < 0 ~ "T-bet down-regulated"
          )) |>
          dplyr::pull(status) |>
          factor()
        limma::plotMA(efit,
          status = status,
          main = paste0("KO vs. CTL MA plot of ", efit$cluster, " DE analysis"),
          values = c("T-bet up-regulated", "T-bet down-regulated", "T-bet"),
          hl.col = c("red", "blue", "blue"),
          hl.pch = c(16, 16, 17),
          hl.cex = c(0.7, 0.7, 1),
          legend = FALSE
        )
        dev.off()
        return(file.path("output", "KO_CTL_MAplot_Bcell.pdf"))
      })()
  ),
  tar_target(KO_CTL_MAplot_all,
    format = "file",
    efits_koctrl |>
      subset(names(efits_koctrl) != "B cell") |>
      (\(efits) {
        pdf(file.path("output", "KO_CTL_MAplot.pdf"), width = 10, height = 16)
        par(mfrow = c(3, 2))
        for (efit in efits) {
          status <- efit$genes |>
            dplyr::left_join(prev_de, by = c("symbol" = "Symbol")) |>
            mutate(status = case_when(
              is.na(logFC) ~ "Other",
              symbol == "Tbx21" ~ "T-bet",
              logFC > 0 ~ "T-bet up-regulated",
              logFC < 0 ~ "T-bet down-regulated"
            )) |>
            dplyr::pull(status) |>
            factor()
          limma::plotMA(efit,
            status = status,
            main = paste0("KO vs. CTL MA plot of ", efit$cluster, " DE analysis"),
            values = c("T-bet up-regulated", "T-bet down-regulated", "T-bet"),
            hl.col = c("red", "blue", "blue"),
            hl.pch = c(16, 16, 17),
            hl.cex = c(0.7, 0.7, 1)
          )
        }
        dev.off()
        return(file.path("output", "KO_CTL_MAplot.pdf"))
      })()
  ),
  tar_target(
    KO_CTL_barcodeplot_all,
    format = "file",
    efits_koctrl |>
      subset(names(efits_koctrl) != "B cell") |>
      (\(efits) {
        pdf(file.path("output", "KO_CTL_barcodeplot.pdf"), width = 14, height = 16)
        par(mfrow = c(3, 2))
        for (efit in efits) {
          status <- efit$genes |>
            dplyr::left_join(prev_de, by = c("symbol" = "Symbol")) |>
            mutate(
              status = case_when(
                is.na(logFC) ~ "Other",
                logFC > 0 ~ "T-bet up-regulated",
                logFC < 0 ~ "T-bet down-regulated"
              )
            ) |>
            dplyr::pull(status) |>
            factor(levels = c("T-bet up-regulated", "T-bet down-regulated", "Other"))
          barcodeplot(efit$t[, "Kovsctrl"],
            index = status == "T-bet up-regulated", index2 = status == "T-bet down-regulated",
            main = paste0("Barcode plot of T-bet signature genes in KO vs. CTL ", efit$cluster, " DE analysis"),
          )
        }
        dev.off()
        return(file.path("output", "KO_CTL_barcodeplot.pdf"))
      })()
  ),
  tar_target(
    KO_CTL_barcodeplot_Bcell,
    format = "file",
    efits_koctrl[["B cell"]] |>
      (\(efit) {
        pdf(file.path("output", "KO_CTL_barcodeplot_Bcell.pdf"), width = 7, height = 5)
        status <- efit$genes |>
          dplyr::left_join(prev_de, by = c("symbol" = "Symbol")) |>
          mutate(
            status = case_when(
              is.na(logFC) ~ "Other",
              logFC > 0 ~ "T-bet up-regulated",
              logFC < 0 ~ "T-bet down-regulated"
            )
          ) |>
          dplyr::pull(status) |>
          factor(levels = c("T-bet up-regulated", "T-bet down-regulated", "Other"))
        barcodeplot(efit$t[, "Kovsctrl"],
          index = status == "T-bet up-regulated", index2 = status == "T-bet down-regulated",
          main = paste0("Barcode plot of T-bet signature genes in KO vs. CTL ", efit$cluster, " DE analysis"),
        )
        dev.off()
        return(file.path("output", "KO_CTL_barcodeplot_Bcell.pdf"))
      })()
  ),

  # enrichment analysis
  ## GO
  tar_map(
    values = tibble::tibble(
      ont = c("BP", "MF", "CC")
    ), names = "ont",
    tar_target(go_koctrl,
      sapply(
        efits_koctrl,
        function(efit) {
          limma::topTable(efit, number = Inf, p.value = 0.05) |>
          filter(!symbol %in% prev_de$Symbol) |>
          rownames() |>
          my_enrichGO(key_type = "ensembl_gene_id", ont = ont)
        }
      )
    ),
    tar_target(go_dotplot_koctrl,
      format = "file",
      go_koctrl |>
        (\(gos) {
          pdf(file.path("output", paste0("KO_CTL_GO_", ont, "_dotplot.pdf")), width = 10, height = 16)
          for (i in names(gos)) {
            if (nrow(head(gos[[i]])) == 0) {
              next
            }
            clusterProfiler::dotplot(gos[[i]],
              title = paste0("KO vs. CTL ", i, " GO ", ont, " enrichment")
            ) |>
              plot()
          }
          dev.off()
          return(file.path("output", paste0("KO_CTL_GO_", ont, "_dotplot.pdf")))
        })()
    ),
    tar_target(go_goplot_koctrl,
      format = "file",
      go_koctrl |>
        (\(gos) {
          pdf(file.path("output", paste0("KO_CTL_GO_", ont, "_goplot.pdf")), width = 25, height = 16)
          par(mfrow = c(3, 2))
          for (i in names(gos)) {
            if (nrow(head(gos[[i]])) == 0) {
              next
            }
            plot(clusterProfiler::goplot(gos[[i]]) + 
            ggtitle(paste0("KO vs. CTL ", i, " GO ", ont, " enrichment")))
          }
          dev.off()
          return(file.path("output", paste0("KO_CTL_GO_", ont, "_goplot.pdf")))
        })()
    ),
    tar_target(go_koctrl_results,
      format = "file",
      {
        go_koctrl |>
          lapply(function(x) {
            x@result |>
              filter(p.adjust < 0.05)
          }) |> 
          bind_rows(.id = "cluster") |>
          write_csv(file.path("output", paste0("KO_CTL_GO_", ont, "_results.csv")))
          file.path("output", paste0("KO_CTL_GO_", ont, "_results.csv"))
      }
    )
  ),

  # KEGG
  tar_target(
    kegg_koctrl,
    sapply(
      efits_koctrl,
      function(efit) {
        limma::topTable(efit, number = Inf, p.value = 0.05) |>
          filter(!symbol %in% prev_de$Symbol) |>
          rownames() |>
          my_enrichKEGG(key_type = "ensembl_gene_id")
      }
    )
  ),
  tar_target(
    kegg_dotplot_koctrl,
    format = "file",
    kegg_koctrl |>
      (\(keggs) {
        pdf(file.path("output", "KO_CTL_KEGG_dotplot.pdf"), width = 10, height = 16)
        for (i in names(keggs)) {
          clusterProfiler::dotplot(keggs[[i]],
            title = paste0("KO vs. CTL ", i, " KEGG enrichment")
          ) |>
            plot()
        }
        dev.off()
        return(file.path("output", "KO_CTL_KEGG_dotplot.pdf"))
      })()
  ),
  tar_target(
    kegg_koctrl_results,
    format = "file",
    {
      kegg_koctrl |>
        lapply(function(x) {
          x@result |>
            filter(p.adjust < 0.05)
        }) |>
        bind_rows(.id = "cluster") |>
        write_csv(file.path("output", "KO_CTL_KEGG_results.csv")
        )
      file.path("output", "KO_CTL_KEGG_results.csv")
    }
  )
)
