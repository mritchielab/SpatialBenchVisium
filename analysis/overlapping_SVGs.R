# SVG overlaps
# Author: Mei Du

library(here)
library(SpatialExperiment)
library(tidyverse)

oct <- readRDS(here("analysis","output","RDS","FF_SVGs_sig.RDS"))
ffpe <- readRDS(here("analysis","output","RDS","FFPE_SVGs_sig.RDS"))
oct_ca <- readRDS(here("analysis","output","RDS","FF_CA_SVGs_sig.RDS"))
ffpe_ca <- readRDS(here("analysis","output","RDS","FFPE_CA_SVGs_sig.RDS"))

#---------------------------Supplementary Figure S3b-----------------------------
## 2-way for FFPE and FFPE CA
library(ggvenn)
venn_FFPE <- list(FFPE=ffpe$gene_id, FFPE_CA=ffpe_ca$gene_id) %>%
  ggvenn(
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 6, text_size=7
  )

## 2-way for OCT and OCT CA
venn_OCT <- list(OCT=oct$gene_id, OCT_CA=oct_ca$gene_id) %>%
  ggvenn(
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 6, text_size=7
  )

pdf(here("analysis","output","svg_overlap_venn.pdf"), height=8, width=6)
cowplot::plot_grid(venn_FFPE,venn_OCT, nrow = 2)
dev.off()

#----------------------------------Figure 3d------------------------------------
venn_FFPE <- list(FFPE=ffpe$gene_id[1:1000], FFPE_CA=ffpe_ca$gene_id[1:1000]) %>%
  ggvenn(
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 6, text_size=7
  )

## c. 2-way for OCT and OCT CA
venn_OCT <- list(OCT=oct$gene_id[1:1000], OCT_CA=oct_ca$gene_id[1:1000]) %>%
  ggvenn(
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 6, text_size=7
  )

pdf(here("analysis","output","top1000_svg_overlap_venn.pdf"), height=8, width=6)
cowplot::plot_grid(venn_FFPE,venn_OCT, nrow = 2)
dev.off()

#----------------------------Supplementary Figure S3a----------------------------
library(ComplexUpset)
list_input <- list(OCT=oct$gene_id, OCT_CA = oct_ca$gene_id, FFPE=ffpe$gene_id, FFPE_CA=ffpe_ca$gene_id)
colSums(UpSetR::fromList(list_input)) # set sizes
pdf(here("analysis","output","svg_overlap_upset.pdf"), height=6, width=10)
ComplexUpset::upset(UpSetR::fromList(list_input), intersect=c("OCT","OCT_CA","FFPE","FFPE_CA"),
               sort_sets=FALSE,
               base_annotations = list(
                 'Intersection size'=(
                   intersection_size()
                   + theme(plot.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.title = element_text(size=17),
                           axis.text = element_text(size=17)))),
               themes=upset_modify_themes(
                 list('intersections_matrix'=theme(text=element_text(size=17), axis.title.x = element_blank()),
                      'overall_sizes'=theme(axis.text.x=element_text(size=15), axis.title.x = element_text(size=17)))))
dev.off()

#----------------------------------Figure 3c------------------------------------
list_input <- list(OCT=oct$gene_id[1:1000], OCT_CA = oct_ca$gene_id[1:1000], FFPE=ffpe$gene_id[1:1000], FFPE_CA=ffpe_ca$gene_id[1:1000])
# pdf(here("analysis","output","top1000_svg_overlap_upset.pdf"), height=6, width=11)
ComplexUpset::upset(UpSetR::fromList(list_input), intersect=c("OCT","OCT_CA","FFPE","FFPE_CA"),
               sort_sets=FALSE,
               base_annotations = list(
                 'Intersection size'=(
                   intersection_size()
                   + theme(plot.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.title = element_text(size=17),
                           axis.text = element_text(size=17)))),
               themes=upset_modify_themes(
                 list('intersections_matrix'=theme(text=element_text(size=17), axis.title.x = element_blank()),
                      'overall_sizes'=theme(axis.text.x=element_text(size=15), axis.title.x = element_text(size=17)))))
# dev.off()

#----------------------------Supplementary Figure S3c----------------------------
source(here("analysis","FFPE_CA_multi-sample.Rmd"))
upset(fromList(top),
      order.by = "freq",
      keep.order=TRUE,
      nsets=5,
      set_size.show=FALSE,
      text.scale=c(1.5, 1.5, 1.5, 1.5, 1.5, 2))

#----------------------------Supplementary Figure S3d----------------------------
dotplot(go) + scale_colour_continuous(labels=scales::scientific_format(digits=1),low="royalblue3",high="orange")


