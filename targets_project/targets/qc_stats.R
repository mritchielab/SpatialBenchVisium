list(
  # tar_target(
  #   all_samples,
  #   load_samples_tb(sample_name_func = function(x) {
  #     gsub("^[^_]+_", "", x$folder)
  #   })
  # ),
  tar_target(all_samples, qs::qread("data/all_samples.qs")),
  tar_target(
    umi_counts,
    all_samples |>
      (\(x){
        return(data.frame(
          umi_all = x$raw_spe |> sapply(\(spe) sum(counts(spe))),
          umi_in_tissue = x$raw_spe |> sapply(\(spe) sum(counts(spe[, spe$in_tissue]))),
          sample_names = x$sample
        ))
      })()
  ),
  tar_target(
    polyA_gene_counts,
    all_samples |>
      filter(experiment == "OCT Experiment 2") |>
      pull("raw_spe") |>
      sapply(\(x) x[, x$in_tissue]) |>
      # genes with a count of 3 or more in at least 10% of spots under tissue
      sapply(\(x) scuttle::addPerFeatureQC(x, threshold = 3)) |>
      sapply(\(x) rownames(x[rowData(x)$detected > 10, ])) |>
      # merge unique names
      (\(x) Reduce(union, x))() |>
      length()
  ),
  tar_target(
    prob_based_gene_counts,
    all_samples |>
      filter(experiment %in% c("CA Experiment 4", "FFPE Experiment 1", "FFPE Experiment 3")) |>
      pull("raw_spe") |>
      sapply(\(x) x[, x$in_tissue]) |>
      sapply(\(x) scuttle::addPerFeatureQC(x, threshold = 3)) |>
      sapply(\(x) rownames(x[rowData(x)$detected > 10, ])) |>
      (\(x) Reduce(union, x))() |>
      length()
  ),
  tar_target(
    prob_based_gene_counts_by_experiment,
    all_samples |>
      filter(experiment %in% c("CA Experiment 4", "FFPE Experiment 1", "FFPE Experiment 3")) |>
      mutate(new_group = case_when(
        experiment %in% c("FFPE Experiment 1", "FFPE Experiment 3") ~ "FFPE",
        experiment == "CytAssist" & group == "WT OCT" ~ "OCT CA",
        experiment == "CytAssist" & group == "WT FFPE" ~ "FFPE CA"
      )) |>
      (\(samples) {
        gene_counts <- c()
        for (i in c("FFPE", "OCT CA", "FFPE CA")) {
          gene_counts[i] <- samples |>
            filter(new_group == i) |>
            pull("raw_spe") |>
            sapply(\(x) x[, x$in_tissue]) |>
            sapply(\(x) scuttle::addPerFeatureQC(x, threshold = 3)) |>
            sapply(\(x) rownames(x[rowData(x)$detected > 10, ])) |>
            (\(x) Reduce(union, x))() |>
            length()
        }
        gene_counts
      })()
  ),
  # umis.pdf
  tar_target(
    umis_pdf,
    format = "file",
    all_samples |>
      (\(all_samples) {
        all_samples$p_umi <- names(all_samples$raw_spe) |>
          sapply(
            \(x) {
              spe <- all_samples$raw_spe[[x]][, all_samples$raw_spe[[x]]$in_tissue]
              if (ncol(spe) < 1000) {
                spot_size <- 1
              } else if (ncol(spe) < 2000) {
                spot_size <- 0.6
              } else {
                spot_size <- 0.4
              }
              tibble(
                counts = colSums(counts(spe)),
                "x_coor" = spatialCoords(spe)[, 1],
                "y_coor" = -spatialCoords(spe)[, 2]
              ) |>
                ggplot(aes(x = x_coor, y = y_coor, col = counts)) +
                geom_point(size = spot_size) +
                coord_fixed() +
                theme_minimal() +
                scale_colour_gradientn(
                  colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
                ) +
                xlab("") +
                ylab("") +
                theme(
                  plot.margin = unit(c(1, 0, 0, 0), "cm"),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank()
                ) +
                ggtitle(gsub("^.*_", "", gsub("(A)|(B)$", "", x)))
            },
            simplify = F
          )

        umi_grid <- plot_grid(
          plot_grid(
            nrow = 1, hjust = 0,
            plotlist = all_samples |>
              filter(experiment == "FFPE Experiment 1") |>
              pull(p_umi)
          ),
          plot_grid(
            nrow = 1, hjust = 0,
            plotlist = all_samples |>
              filter(experiment == "OCT Experiment 2" & group != "WT OCT") |>
              pull(p_umi)
          ),
          plot_grid(
            nrow = 1, hjust = 0,
            plotlist = all_samples |>
              filter(experiment == "OCT Experiment 2" & group == "WT OCT") |>
              pull(p_umi)
          ),
          plot_grid(
            nrow = 1, hjust = 0,
            plotlist = all_samples |>
              filter(experiment == "FFPE Experiment 3") |>
              pull(p_umi)
          ),
          plot_grid(
            plot_grid(
              nrow = 1, hjust = 0, vjust = 6,
              plotlist = all_samples$p_umi[c("CytAssist_OCT_709", "CytAssist_OCT_713")]
            ),
            plot_grid(
              nrow = 1, hjust = 0, vjust = 6,
              plotlist = all_samples$p_umi[c("CytAssist_FFPE_709", "CytAssist_FFPE_713")]
            ),
            nrow = 1, vjust = 2.5,
            labels = c("OCT CA", "FFPE CA")
          ),
          ncol = 1, vjust = 1, hjust = 0,
          labels = c("FFPE manual (earlier)", "OCT manual (KO vs CTRL)", "OCT manual (WT)", "FFPE manual (later)", "CA")
        )

        cowplot::save_plot(
          plot = umi_grid,
          base_height = 14, base_width = 13, filename = file.path("output", "umi.pdf"), dpi = 600
        )
        return(file.path("output", "umi.pdf"))
      })()
  ),

  # MDS
  tar_target(samples_mds, format = 'file',
    {
      mei_path <- "/vast/scratch/users/wang.ch/tmp/SpatialBench"
      svgs <- Reduce(intersect, 
        list(
          readRDS(file.path(mei_path, "output", "FF_SVGs_sig.RDS"))$gene_name[1:1000],
          readRDS(file.path(mei_path, "output", "FFPE_SVGs_sig.RDS"))$gene_name[1:1000],
          readRDS(file.path(mei_path, "output", "FF_CA_SVGs_sig.RDS"))$gene_name[1:1000],
          readRDS(file.path(mei_path, "output", "FFPE_CA_SVGs_sig.RDS"))$gene_name[1:1000]
        )
      )
      spes <- scuttleFilter(filter(all_samples, experiment != 'Pilot study'))$spe
      spes <- lapply(spes, \(x) {rownames(x) <- rowData(x)$symbol; x[svgs, ]})
      bulkCounts <- do.call(cbind, lapply(spes, \(x) rowSums(counts(x))))
      mds <- edgeR::DGEList(
          counts = bulkCounts,
          samples = colnames(bulkCounts),
      ) |>
        edgeR::calcNormFactors() |>
        limma::plotMDS(plot = FALSE)
      tb <- tibble(
        x = mds$x,
        y = mds$y,
        experiment = dplyr::filter(all_samples, experiment != 'Pilot study')$experiment,
        sample = colnames(bulkCounts),
        group = dplyr::filter(all_samples, experiment != 'Pilot study')$group
      ) |>
        mutate(experiment = case_when(
          experiment == "FFPE Experiment 1" ~ "FFPE manual (earlier)",
          experiment == "FFPE Experiment 3" ~ "FFPE manual (later)",
          experiment == "CA Experiment 4" ~ "CytAssist",
          experiment == "OCT Experiment 2" & group == "WT OCT" ~ "OCT manual (WT)",
          experiment == "OCT Experiment 2" & group != "WT OCT" ~ "OCT manual (KOvsCTRL)"
        ), protocol = case_when(
          grepl("FFPE manual", experiment) ~ "FFPE manual",
          grepl("OCT manual", experiment) ~ "OCT manual",
          grepl("CytAssist_OCT", sample) ~ "OCT CA",
          grepl("CytAssist_FFPE", sample) ~ "FFPE CA"
        )) |>
        mutate(experiment = factor(experiment, levels = c("FFPE manual (earlier)", "OCT manual (KOvsCTRL)", "OCT manual (WT)","FFPE manual (later)", "CytAssist")),
               protocol = factor(protocol, levels = c("OCT manual","FFPE manual","OCT CA", "FFPE CA")),
          sample = stringr::str_extract(sample, "[0-9]{3}")
        )

      p <- ggpubr::ggscatter(tb, x="x", y="y",
        color = "experiment", palette = "jama",
        shape = "protocol",  label="sample", repel=TRUE,
        font.label=c(5,"plain"),
        xlab = sprintf("MDS1 (%.2f%%)", mds$var.explained[1]*100),
        ylab = sprintf("MDS2 (%.2f%%)", mds$var.explained[2]*100)) +
        theme(legend.position= c(0.68,0.95),
          legend.justification = c("left","top"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box = "horizontal",legend.text = element_text(size=7),
          legend.title=element_text(size=7),
          axis.title = element_text(face="bold", size=5),
          axis.text=element_text(size=5))
      
      ggsave(file.path("output", "mds_samples.pdf"), p, width = 8.5, height = 6.5)
      return(file.path("output", "mds_samples.pdf"))
    }
  )
)
