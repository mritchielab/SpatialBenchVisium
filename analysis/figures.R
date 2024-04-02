# SpatialBench paper figures
# Author: Mei Du

library(here)
library(SpatialExperiment)
library(Seurat)
library(cowplot)
library(ggspavis)
library(RColorBrewer)
library(scater)
library(scran)
library(pheatmap)
library(tidyverse)
library(Mus.musculus)
library(biomaRt)

dir <- "/stornext/Projects/score/Analyses/G000218_spatial_benchmarking_study/extdata/SpaceRanger"

list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}
files <- list.dirs.depth.n(dir, 2)

#--------------fig 1 data quality----------------
# c) - UMI counts
sample <- c(460,462,463,544,545,708,709,173,174,167,168,708,709,710,713,709,713,709,713)
run <- c(rep("V11L12-124",3),
         rep("V11D13-369",4),
         rep("V12F28-091",4),
         rep("V12M14-014",4),
         rep("V42L05-065",2),
         rep("V42L05-390",2))
protocol <- c(rep("FFPE",3),
           rep("OCT",4),
           rep("OCT",4),
           rep("FFPE",4),
           rep("OCT CA",2),
           rep("FFPE CA",2))
experiment <- c(rep("FFPE (1)", 3),
                rep("OCT (2)",4),
                rep("KOvsCTRL (2)",4),
                rep("FFPE (3)", 4),
                rep("CytAssist (4)",4))

metadata <- as.data.frame(cbind(sample, run, protocol, experiment))

get_spe_counts <- function(i) {
  spe <- read10xVisium(files[str_detect(files,paste0(metadata[i,"run"],"/.*",metadata[i,"sample"]))], type="sparse", data="raw", images="lowres", load=TRUE)
  umi <- data.frame(sum=colSums(counts(spe)),sample=metadata[i,"sample"], run=metadata[i,"run"])
  umi <- rownames_to_column(umi, var="barcode" )
}

umi_counts <- data.frame()
for (i in 1:length(sample)) {
  umi_counts <- rbind(umi_counts, get_spe_counts(i))
}
data <- right_join(metadata,umi_counts)
data$experiment <- factor(data$experiment,levels=c("FFPE (1)", "KOvsCTRL (2)", "OCT (2)","FFPE (3)", "CytAssist (4)"))
data$protocol <- factor(data$protocol, levels = c("OCT","FFPE","OCT CA", "FFPE CA"))

# saveRDS(data, here("analysis","output","rds","umi_counts.RDS"))

library(ggplot2)

# facet by protocol, colour by experiment
a <- data %>%
  ggplot(aes(fill=experiment, y=sum, x=sample)) +
  ggsci::scale_fill_jama() +
  geom_violin() +
  scale_y_log10(labels=scales::comma, breaks = c(0, 1000, 10000, 50000,100000)) +
  xlab("sample") + ylab ("UMI counts") + facet_wrap(vars(protocol), nrow=2, ncol=2, scales="free_x") +
  theme_light() + theme(strip.background = element_rect(fill = "#DCEFFE", color = "#DCEFFE"),
                          strip.text = element_text(color = "#074D6C", face = "bold", size = 15),
                        axis.title = element_text(face="bold",size=16),
                        axis.text = element_text(size=12))
# extract a legend that is laid out horizontally
legend <- get_legend(
  a +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.text = element_text(size=18), legend.title=element_text(size=18))
)


# d) violin plot of detected genes per spot
get_gene_counts <- function(i) {
  spe <- read10xVisium(files[str_detect(files,paste0(metadata[i,"run"],"/.*",metadata[i,"sample"]))], type="sparse", data="raw", images="lowres", load=TRUE)
  spe <- addPerCellQC(spe)
  genes <- data.frame(barcode=rownames(colData(spe)),sum=colData(spe)$detected,sample=metadata[i,"sample"], run=metadata[i,"run"])
  return(genes)
  genes <- column_to_rownames(genes, var="barcode")
}

gene_counts <- data.frame()
for (i in 1:length(sample)) {
  gene_counts <- rbind(gene_counts, get_gene_counts(i))
}
gene_data <- right_join(metadata,gene_counts)
gene_data$experiment <- factor(gene_data$experiment,levels=c("FFPE (1)", "KOvsCTRL (2)", "OCT (2)","FFPE (3)", "CytAssist (4)"))
gene_data$protocol <- factor(gene_data$protocol, levels = c("OCT","FFPE","OCT CA", "FFPE CA"))

# saveRDS(data, here("analysis","output","rds","gene_counts.RDS"))

library(ggplot2)

# facet by protocol, colour by experiment
b <- gene_data %>%
  ggplot(aes(fill=experiment, y=sum, x=sample)) +
  ggsci::scale_fill_jama() +
  geom_violin() +
  scale_y_continuous(labels=scales::comma, breaks = c(0, 5000, 10000)) +
  xlab("sample") + ylab ("Number of genes per spot") + facet_wrap(vars(protocol), nrow=2, ncol=2, scales="free_x") +
  theme_light() + theme(strip.background = element_rect(fill = "#DCEFFE", color = "#DCEFFE"),
                        strip.text = element_text(color = "#074D6C", face = "bold", size = 15),
                        axis.title = element_text(face="bold",size=16),
                        axis.text = element_text(size=12))

# e) - Mean reads per spot vs fraction of reads in spots under tissue

get_metrics <- function(i) {
  csv <- read.csv(paste0(files[str_detect(files,paste0(metadata[i,"run"],"/.*",metadata[i,"sample"]))],"/outs/metrics_summary.csv"))
  metrics <- data.frame(sample=csv[,"Sample.ID"], mean_reads=csv[,"Mean.Reads.per.Spot"],
                        fraction=csv[,"Fraction.Reads.in.Spots.Under.Tissue"], run = metadata[i,"run"])
}
metrics <- data.frame()
for (i in 1:length(sample)) {
  metrics <- rbind(metrics, get_metrics(i))
}
metrics$sample <- str_extract(metrics$sample, "\\d+[A-Za-z]*$")
metrics$sample <- str_replace(metrics$sample, "[ABC]","")
metrics_full <- right_join(metadata,metrics)
metrics_full$experiment <- factor(metrics_full$experiment,levels=c("FFPE (1)", "KOvsCTRL (2)", "OCT (2)","FFPE (3)", "CytAssist (4)"))
metrics_full$protocol <- factor(metrics_full$protocol, levels = c("OCT","FFPE","OCT CA", "FFPE CA"))

library(ggpubr)
c<-ggscatter(metrics_full,
             x="fraction",
             y="mean_reads",
             color = "experiment", palette = "jama",
             shape = "protocol" , size=5, label="sample", repel=TRUE,
             font.label=c(15,"plain"),
             xlab = "Fraction of reads in spots under tissue",
             ylab = "Mean reads per spot") +theme(
               legend.position = c(0.7,0.95),
               legend.justification = c("left","top"),
               legend.margin = margin(6, 6, 6, 6),
               legend.box = "horizontal",legend.text = element_text(size=17),
               legend.title=element_text(size=17),
               axis.title = element_text(face="bold", size=15),
               axis.text=element_text(size=14)
             ) + guides(color = "none")


# add the legend underneath the row we made earlier. Give it 10% of the height of one plot (via rel_heights).
row <- plot_grid(a+theme(legend.position="none"),b+theme(legend.position="none"), c,labels=c("c","d","e"), nrow=1,scale=c(1,1,0.95),label_size=30)

pdf(here("analysis","output","data_quality_paper.pdf"), height=8,width=22)
cowplot::plot_grid(legend,row,ncol=1, rel_heights = c(.1, 1), scale=c(2,1))
dev.off()


#------------fig 2 poly A vs probe----------------
# b)
source(here("analysis", "EDA_709_FFPE_CA.Rmd"))
# pdf(here("analysis","output","FFPE_CA_709_UMI_violin.pdf"),height=4.5,width=5)
umi_violin
# dev.off()
# pdf(here("analysis","output","FFPE_CA_709_genes_filtered.pdf"),height=5,width=5)
genes_detected
# dev.off()
# pdf(here("analysis","output","FFPE_CA_709_mito_filtered.pdf"),height=5,width=5)
mito_content
# dev.off()

# d)
source(here("analysis","overlapping_genes.R"))
# pdf(here("analysis","output","fig2_de.pdf"), height=6,width=18)
plot_grid(upset, venn_all, rel_widths=c(1.3,0.9), scale=c(1,1),labels=c("a","b"),  label_size=30)
# dev.off()

#------------fig 3 downstream analysis------------
source(here("analysis","FFPE_CA_multi-sample.Rmd"))

# a) HVG and SVG spatial expression
# bottom
# pdf(here("analysis","output","FFPE_CA_top2_SVGs.pdf"),height=2.3, width=7)
p1+p2
# dev.off()

# top - use ix_gene <- c("Car2","Tmcc2") in "Spatial SVG expression"

# b) clustering using single-cell methods
# pdf(here("analysis","output","FFPE_CA_709_cluster_UMAP.pdf"),height=4,width=5)
umap_709
# dev.off()
source(here("analysis", "EDA_713_FFPE_CA.Rmd"))
# pdf(here("analysis","output","FFPE_CA_713_cluster_UMAP.pdf"),height=4,width=5)
umap_713
# dev.off()

# c), d) - see "overlapping_SVGs.R"

# e) - multi-sample clustering
umap2
umap3 # supp fig S6
# e) bottom, T cell cluster
c7
# f) cluster score heatmap
p
# g) deconvolution
deconv_plot
# h) annotated clusters
labelled_plot

#------------supplementary figure S4------------
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15

#------------supplementary figure S5------------
source(here("analysis","qupath_classification.R"))
plot_709 + plot_713
