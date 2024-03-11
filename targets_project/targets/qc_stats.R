list(
  tar_target(
    all_samples,
    load_samples_tb(sample_name_func = function(x) {
      gsub("^[^_]+_", "", x$folder)
    })
  ),
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
          labels = c("FFPE Experiment 1 (WT)", "OCT Experiment 2 (KO & CTRL)", "OCT Experiment 2 (WT)", "FFPE Experiment 3 (WT)", "CA Experiment 4 (WT)")
        )

        cowplot::save_plot(
          plot = umi_grid,
          base_height = 14, base_width = 13, filename = file.path("output", "umi.pdf"), dpi = 600
        )
        return(file.path("output", "umi.pdf"))
      })()
  )
)
