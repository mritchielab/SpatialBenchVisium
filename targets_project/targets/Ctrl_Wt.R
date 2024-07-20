#  WT vs CTRL
list(
  tar_target(spes_wt_ctrl_tb, 
    all_samples |>
      filter(experiment == "OCT Experiment 2" & group != "KO OCT") |> 
      mutate(sample = gsub(".*_", "", folder)) |>
      scuttleFilter(feature_discard_func = function(x) {
        FALSE
      })
  ),
  tar_target(
    spe_wtctrl,
    spes_wt_ctrl_tb |>
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
      (\(x) { # add sex labels
        x$sex <- factor(ifelse(x$sample_id %in% c("173", "174"), "male", "female"))
        return(x)
      })()
  ),
  tar_target(iSCMEB_wtctrl, 
    {
      set.seed(123)
      if (file.exists('output/iSCMEB_wtctrl.qs')) {
        qs::qread('output/iSCMEB_wtctrl.qs')
      } else {
        spes_wt_ctrl_tb$spe |>
          sapply(to_seu) |>
          run_iSC_MEB(customGenelist = filter_nnSVG(run_nnSVG(spes_wt_ctrl_tb$spe)), k = 8)
      }
    }
  ),
  tar_target(
    cluster_wtctrl,
    iSCMEB_wtctrl@resList@idents |>
      unlist() |>
      factor()
  ),
  tar_target(
    plot_umap_wtctrl,
    format = "file",
    {
      p <- iSCMEB_wtctrl |>
        iSC.MEB::CalculateUMAP(reduction = "iSCMEB", n_comp = 2) |>
        iSC.MEB::LowEmbedPlot(item = "cluster", reduction = "UMAP2", point_size = 0.1, point_alpha = 0.8) +
        theme(
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20)
        ) +
        guides(color = guide_legend(override.aes = list(size = 4)))
      ggsave(file.path("output", "WT_CTL_UMAP.pdf"), p, width = 15, height = 12, dpi = 600, units = "cm")
      file.path("output", "WT_CTL_UMAP.pdf")
    }
  ),
  tar_target(
    plot_zone_wtctrl,
    format = "file",
    {
      tb <- spe_wtctrl |>
        spatialCoords() |>
        data.frame() |>
        as_tibble() |>
        mutate(
          sex = spe_wtctrl$sex,
          sample_id = spe_wtctrl$sample_id,
          cluster = zone_wtctrl
        )

      theme <- theme_minimal() +
        theme(
          plot.margin = unit(c(1, 0, 0, 0), "cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)
        )

      p1 <- tb |>
        filter(sex == "male") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.7) +
        coord_fixed() +
        facet_wrap(~sample_id, nrow = 1) +
        theme +
        scale_color_manual(values = mei_palette) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        xlab("") +
        ylab("") +
        ggtitle("Male")
      p2 <- tb |>
        filter(sex == "female") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.5) +
        coord_fixed() +
        facet_wrap(~sample_id, nrow = 1) +
        theme +
        scale_color_manual(values = mei_palette) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        xlab("") +
        ylab("") +
        ggtitle("Female")

      p <- ((p1 + guides(color = "none")) / p2) + plot_layout(guides = "collect")
      ggsave(file.path("output", "WT_CTL_zone_plot.pdf"), p, width = 30, height = 17.5, dpi = 600, units = "cm")
      file.path("output", "WT_CTL_zone_plot.pdf")
    }
  ),
  tar_target(
    plot_cluster_wtctrl,
    format = "file",
    {
      tb <- spe_wtctrl |>
        spatialCoords() |>
        data.frame() |>
        as_tibble() |>
        mutate(
          sex = spe_wtctrl$sex,
          sample_id = spe_wtctrl$sample_id,
          cluster = cluster_wtctrl
        )

      theme <- theme_minimal() +
        theme(
          plot.margin = unit(c(1, 0, 0, 0), "cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)
        )

      p1 <- tb |>
        filter(sex == "male") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.7) +
        coord_fixed() +
        facet_wrap(~sample_id, nrow = 1) +
        theme +
        scale_color_manual(values = iSC.MEB_palette(length(levels(cluster_wtctrl)))) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        xlab("") +
        ylab("") +
        ggtitle("Male")
      p2 <- tb |>
        filter(sex == "female") |>
        ggplot(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, col = cluster)) +
        geom_point(size = 0.5) +
        coord_fixed() +
        facet_wrap(~sample_id, nrow = 1) +
        theme +
        scale_color_manual(values = iSC.MEB_palette(length(levels(cluster_wtctrl)))) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        xlab("") +
        ylab("") +
        ggtitle("Female")

      p <- ((p1 + guides(color = "none")) / p2) + plot_layout(guides = "collect")
      ggsave(file.path("output", "WT_CTL_cluster_plot.pdf"), p, width = 30, height = 17.5, dpi = 600, units = "cm")
      file.path("output", "WT_CTL_cluster_plot.pdf")
    }
  ),

  # marker heatmap
  tar_target(zone_markers_df, list(
    "Macrophage" = c("Cd274", "Marco", "Csf1r", "Adgre1", "Cd209b", "Cd206", "Cd80", "Mac1", "Cd68")[c(1, 3, 4, 5, 7, 9)],
    "B cell" = c("Cd19", "Cd22", "Ighd", "Cd5"),
    "Germinal centre" = c("Cxcr4", "Cd83", "Bcl6", "Rgs13", "Aicda"),
    "Neutrophil" = c("S100a9", "S100a8", "Ngp"),
    "Erythrocyte" = c("Car2", "Car1", "Klf1"),
    "Plasma cell" = c("Cd38", "Cd138", "Xbp1", "Irf4", "Prdm1", "Cd27", "Cd319", "Mum1")[c(1, 3, 4, 5, 6)],
    "T cell" = c("Trac", "Cd3d", "Cd4", "Cd3e", "Cd8a"),
    "Red pulp" = c("Ifitm3", "C1qc", "Hmox1", "Hba-a1", "Klf1"),
    "Marginal zone" = c("Marco", "Lyz2", "Ighd", "Igfbp7", "Igfbp3", "Ly6d"),
    "White pulp" = c("Ighd", "Cd19", "Trac", "Trbc2")
  ) |>
    stack() |>
    setNames(c("gene", "zone")) |>
    merge(rowData(spe_wtctrl)[, c("symbol", "EnsembleID")],
      by.x = "gene",
      by.y = "symbol"
    ) |>
    na.omit()),
  tar_target(mei_palette, RColorBrewer::brewer.pal(6, "Set2")[c(1, 3, 5, 4, 6, 2)] |>
    setNames(c("B cell", "Germinal centre", "Neutrophil", "Erythrocyte", "Plasma cell", "T cell"))),
  tar_target(
    heatmap_wtctrl,
    format = "file",
    {
      p <- logfc_mtx(spe_wtctrl, cluster_wtctrl, zone_markers_df) |>
        ggplot(aes(
          x = factor(cluster, levels = levels(cluster_wtctrl)),
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
      ggsave(file.path("output", "WT_CTL_heatmap.pdf"), p, width = 15, height = 15, dpi = 600, units = "cm")
      file.path("output", "WT_CTL_heatmap.pdf")
    }
  ),

  # https://doi.org/10.1371%2Fjournal.pgen.1005079
  tar_target(
    xie_genes,
    c(
      "Cybb", "Ddx3x", "Kdm6a", "Cfp", "Utp14a", "Firre", "Bgn", "5430427O19Rik",
      "Eif2s3x", "Vsig4", "Xist", "Ftx", "5530601H04Rik", "Pbdc1", "5730416F02Rik",
      "Kdm5c", "Tmsb4x"
    )
  ),

  # manual annotation for cluster
  tar_target(zone_wtctrl, cluster_wtctrl |>
    (\(x) {
      x$zone <- factor(
        case_when(
          x == 1 ~ "Neutrophil",
          x == 2 ~ "Germinal centre",
          x == 3 ~ "Erythrocyte",
          x == 4 ~ "Erythrocyte",
          x == 5 ~ "B cell",
          x == 6 ~ "Erythrocyte",
          x == 7 ~ "Plasma cell",
          x == 8 ~ "T cell"
        )
      )
    })()),

  # pseudo-bulk DE
  tar_target(
    puseudo_bulk_counts_wtctrl,
    sapply(levels(zone_wtctrl),
      function(x) {
        aggregate_sample_bulk_counts(spe_wtctrl[, zone_wtctrl == x], set.names = "EnsembleID")
      },
      simplify = F
    )
  ),
  tar_target(
    vooms_wtctrl,
    sapply(puseudo_bulk_counts_wtctrl,
      function(x) {
        dgelist <- DGEList(
          counts = x,
          samples = colnames(x),
          gene = rowData(spe_wtctrl)[
            match(rownames(x), rowData(spe_wtctrl)$EnsembleID),
            c("symbol", "chromosome_name")
          ],
          group = factor(ifelse(colnames(x) %in% c("173", "174"), "male", "female"))
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
    efits_wtctrl,
    sapply(names(vooms_wtctrl),
      function(x) {
        efit <- lmFit(vooms_wtctrl[[x]], vooms_wtctrl[[x]]$design) |>
          contrasts.fit(contrasts = makeContrasts(
            MalevsFemale = "groupmale - groupfemale",
            levels = vooms_wtctrl[[x]]$design
          )) |>
          eBayes()
        efit$voom.line <- vooms_wtctrl[[x]]$voom.line
        efit$voom.xy <- vooms_wtctrl[[x]]$voom.xy
        efit$cluster <- x
        return(efit)
      },
      simplify = F
    )
  ),
  tar_target(
    # male vs female, Y as +1, XiE as -1, expect gene set up-regulated
    roasts_wtctrl,
    sapply(levels(zone_wtctrl),
      function(x) {
        roast(vooms_wtctrl[[x]],
          design = efits_wtctrl[[x]]$design,
          index = efits_wtctrl[[x]]$gene %>%
            mutate(idx = seq.int(nrow(.))) %>%
            filter(chromosome_name == "Y" | symbol %in% xie_genes) %>%
            pull(idx),
          gene.weights = efits_wtctrl[[x]]$gene %>%
            mutate(weight = ifelse(chromosome_name == "Y", 1, -1)) %>%
            pull(weight),
          contrast = efits_wtctrl[[x]]$contrasts,
          nrot = 99999
        )
      },
      simplify = F
    )
  ),
  tar_target(
    roast_table_wtctrl,
    sapply(roasts_wtctrl, function(x) {
      x$p.value["Up", ]
    }) |>
      t() |>
      cbind(Number.of.genes = sapply(roasts_wtctrl, function(x) {
        x$ngenes.in.set
      }))
  ),

  # plots for all clusters
  tar_target(WT_CTL_sex_MAplot_all,
    format = "file",
    efits_wtctrl |>
      subset(names(efits_wtctrl) != "B cell") |>
      (\(efits) {
        pdf(file.path("output", "WT_CTL_sex_MAplot.pdf"), width = 10, height = 16)
        par(mfrow = c(3, 2))
        for (efit in efits) {
          chr <- efit$gene$chromosome_name
          chr <- ifelse(efit$gene$symbol %in% xie_genes, "XiE", chr)
          chr[!chr %in% c("XiE", "Y")] <- "Other"
          chr <- factor(chr,
            levels = c("XiE", "Y", "Other"),
            labels = c("XiE gene", "chrY gene", "Other")
          )
          limma::plotMA(efit,
            status = chr,
            coef = "MalevsFemale",
            main = paste0("Male vs. Female\nMA plot of ", efit$cluster, " DE analysis"),
            hl.col = c("blue", "red", "black")
          )
        }
        dev.off()
        return(file.path("output", "WT_CTL_sex_MAplot.pdf"))
      })()
  ),
  tar_target(
    WT_CTL_sex_barcodeplot_all,
    format = "file",
    efits_wtctrl |>
      subset(names(efits_wtctrl) != "B cell") |>
      (\(efits) {
        pdf(file.path("output", "WT_CTL_sex_barcodeplot.pdf"), width = 14, height = 16)
        par(mfrow = c(3, 2))
        for (efit in efits) {
          chr <- efit$gene$chromosome_name
          chr <- ifelse(efit$gene$symbol %in% xie_genes, "XiE", chr)
          chr[!chr %in% c("XiE", "Y")] <- "Other"
          chr <- factor(chr,
            levels = c("XiE", "Y", "Other")
          )
          barcodeplot(efit$t[, "MalevsFemale"],
            index = chr == "Y", index2 = chr == "XiE",
            main = paste0("Barcode plot of sex-specific genes in Male vs. Female ", efit$cluster, " DE analysis")
          )
        }
        dev.off()
        return(file.path("output", "WT_CTL_sex_barcodeplot.pdf"))
      })()
  ),

  # individual plots
  tar_target(WT_CTL_sex_MAplot,
    format = "file",
    (\(efit) {
      file_name <- paste0(file.path("output", "WT_CTL_sex_MAplot_"), gsub(" ", "_", efit$cluster), ".pdf")
      chr <- efit$gene$chromosome_name
      chr <- ifelse(efit$gene$symbol %in% xie_genes, "XiE", chr)
      chr[!chr %in% c("XiE", "Y")] <- "Other"
      chr <- factor(chr,
        levels = c("XiE", "Y", "Other"),
        labels = c("XiE gene", "chrY gene", "Other")
      )
      pdf(file_name, width = 8, height = 8)
      limma::plotMA(efit,
        status = chr,
        coef = "MalevsFemale",
        main = paste0("Male vs. Female\nMA plot of ", efit$cluster, " DE analysis"),
        values = c("XiE gene", "chrY gene"),
        hl.col = c("blue", "red")
      )
      dev.off()
      return(file_name)
    })(efits_wtctrl[[1]]),
    pattern = map(efits_wtctrl),
    iteration = "list"
  ),
  tar_target(WT_CTL_sex_barcodeplot,
    format = "file",
    (\(efit) {
      file_name <- paste0(file.path("output", "WT_CTL_sex_barcodeplot_"), gsub(" ", "_", efit$cluster), ".pdf")
      chr <- efit$gene$chromosome_name
      chr <- ifelse(efit$gene$symbol %in% xie_genes, "XiE", chr)
      chr[!chr %in% c("XiE", "Y")] <- "Other"
      chr <- factor(chr,
        levels = c("XiE", "Y", "Other")
      )
      pdf(file_name, width = 9, height = 9)
      barcodeplot(efit$t[, "MalevsFemale"],
        index = chr == "Y", index2 = chr == "XiE",
        main = paste0("Barcode plot of sex-specific genes in Male vs. Female ", efit$cluster, " DE analysis")
      )
      dev.off()
      return(file_name)
    })(efits_wtctrl[[1]]),
    pattern = map(efits_wtctrl),
    iteration = "list"
  )
)
