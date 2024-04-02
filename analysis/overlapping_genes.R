# Compiling a “full set of detected genes” for all WT samples
# Author: Mei Du

library(here)
library(SpatialExperiment)
library(tidyverse)

# FFPE
spe_713_FFPE <- readRDS(here("analysis","output","RDS","spe_713_FFPE_norm.RDS"))
spe_709_FFPE <- readRDS(here("analysis","output","RDS","spe_709_FFPE_norm.RDS"))
spe_710 <- readRDS(here("analysis","output","RDS","spe_710_norm.RDS"))
spe_708_FFPE <- readRDS(here("analysis","output","RDS","spe_708_FFPE_norm.RDS"))
spe_460 <- readRDS(here("analysis","output","RDS","spe_460_norm.RDS"))
spe_462 <- readRDS(here("analysis","output","RDS","spe_462_norm.RDS"))
spe_463 <- readRDS(here("analysis","output","RDS","spe_463_norm.RDS"))

# FFPE CA
spe_713_FFPE_CA <- readRDS(here("analysis","output","RDS","spe_713_FFPE_CA_norm.RDS"))
spe_709_FFPE_CA <- readRDS(here("analysis","output","RDS","spe_709_FFPE_CA_norm.RDS"))

# FF
spe_709_FF <- readRDS(here("analysis","output","RDS","spe_709_FF_norm.RDS"))
spe_708_FF <- readRDS(here("analysis","output","RDS","spe_708_FF_norm.RDS"))
spe_545 <- readRDS(here("analysis","output","RDS","spe_545_norm.RDS"))
spe_544 <- readRDS(here("analysis","output","RDS","spe_544_norm.RDS"))
# spe_543 <- readRDS(here("analysis","output","RDS","spe_543_norm.RDS"))
spe_167 <- readRDS(here("analysis","output","RDS","spe_167a_norm.RDS"))
spe_168 <- readRDS(here("analysis","output","RDS","spe_168a_norm.RDS"))
spe_173 <- readRDS(here("analysis","output","RDS","spe_173b_norm.RDS"))
spe_174 <- readRDS(here("analysis","output","RDS","spe_174a_norm.RDS"))

# FF CA
spe_709_FF_CA <- readRDS(here("analysis","output","RDS","spe_709_FF_CA_norm.RDS"))
spe_713_FF_CA <- readRDS(here("analysis","output","RDS","spe_713_FF_CA_norm.RDS"))

#------------1. For each sample, get a list of its detected genes.--------------
# filter out any background noise
# filter_genes_pcspots = 10 for figure 2d-e
# filter_genes_pcspots = 1 for supplementary figure S2
filter <- function(spe) {
  nnSVG::filter_genes(
    spe,
    filter_genes_ncounts = 3,
    filter_genes_pcspots = 10,
    filter_mito = TRUE
  )
}

colnames(rowData(spe_713_FFPE)) <- c("gene_name","Chr")
colnames(rowData(spe_709_FFPE)) <- c("gene_name","Chr")
colnames(rowData(spe_710)) <- c("gene_name","Chr")
colnames(rowData(spe_708_FFPE)) <- c("gene_name","Chr")
colnames(rowData(spe_460)) <- c("gene_name","Chr")
colnames(rowData(spe_462)) <- c("gene_name","Chr")
colnames(rowData(spe_463)) <- c("gene_name","Chr")
spe_713_FFPE_filtered <- filter(spe_713_FFPE) # 5% (n=81), removed 26299/32576 genes
spe_709_FFPE_filtered <- filter(spe_709_FFPE) # 5% (n=92), removed 26294/32576 genes
spe_710_filtered <- filter(spe_710) # 5% (n=94), removed 26601/32576 genes
spe_708_FFPE_filtered <- filter(spe_708_FFPE) # 5% (n=89), removed 26505/32576 genes
spe_460_filtered <- filter(spe_460) # 5% (n=30), removed 25032/32576 genes
spe_462_filtered <- filter(spe_462) # 5% (n=32), removed 25034/32576 genes
spe_463_filtered <- filter(spe_463) # 5% (n=39), removed 25694/32576 genes

colnames(rowData(spe_713_FFPE_CA)) <- c("gene_name","Chr")
colnames(rowData(spe_709_FFPE_CA)) <- c("gene_name","Chr")
spe_713_FFPE_CA_filtered <- filter(spe_713_FFPE_CA) # 5% (n=85), removed 28337/32576 genes
spe_709_FFPE_CA_filtered <- filter(spe_709_FFPE_CA) # 5% (n=97), removed 24782/32576 genes

colnames(rowData(spe_709_FF)) <- c("gene_name","Chr")
colnames(rowData(spe_708_FF)) <- c("gene_name","Chr")
colnames(rowData(spe_545)) <- c("gene_name","Chr")
colnames(rowData(spe_544)) <- c("gene_name","Chr")
colnames(rowData(spe_167a)) <- c("gene_name","Chr")
colnames(rowData(spe_168a)) <- c("gene_name","Chr")
colnames(rowData(spe_173)) <- c("gene_name","Chr")
colnames(rowData(spe_174)) <- c("gene_name","Chr")
spe_709_FF_filtered <- filter(spe_709_FF) # 5% (n=90), removed 30257/32272 genes
spe_708_FF_filtered <- filter(spe_708_FF) # 5% (n=75), removed 31285/32272 genes
spe_545_filtered <- filter(spe_545) # 5% (n=83), removed 30315/32272 genes
spe_544_filtered <- filter(spe_544) # 5% (n=135), removed 29965/32272 genes
spe_167_filtered <- filter(spe_167)
spe_168_filtered <- filter(spe_168)
spe_173_filtered <- filter(spe_173)
spe_174_filtered <- filter(spe_174)

colnames(rowData(spe_709_FF_CA)) <- c("gene_name","Chr")
colnames(rowData(spe_713_FF_CA)) <- c("gene_name","Chr")
spe_709_FF_CA_filtered <- filter(spe_709_FF_CA) # 5% (n=116), removed 25567/32272 genes
spe_713_FF_CA_filtered <- filter(spe_713_FF_CA) # 5% (n=115), removed 25968/32272 genes

#------------2. For each group, combine lists from its samples to get a list of its detected genes. (Groups are: FFPE, FFPE CA, OCT, OCT CA)--------------
FFPE_genes <- unique(c(rownames(spe_713_FFPE_filtered), rownames(spe_709_FFPE_filtered), rownames(spe_710_filtered), rownames(spe_708_FFPE_filtered),
                       rownames(spe_460_filtered), rownames(spe_462_filtered), rownames(spe_463_filtered)))
FFPE_CA_genes <- unique(c(rownames(spe_713_FFPE_CA_filtered), rownames(spe_709_FFPE_CA_filtered)))
OCT_genes <- unique(c(rownames(spe_709_FF_filtered), rownames(spe_708_FF_filtered), rownames(spe_545_filtered), rownames(spe_544_filtered),
                      rownames(spe_167_filtered),rownames(spe_168_filtered),rownames(spe_173_filtered),rownames(spe_174_filtered)))
OCT_CA_genes <- unique(c(rownames(spe_709_FF_CA_filtered), rownames(spe_713_FF_CA_filtered)))
genes_grouped <- list(FFPE=FFPE_genes, FFPE_CA = FFPE_CA_genes, OCT=OCT_genes, OCT_CA=OCT_CA_genes)

#------------3. For the experiment, combine lists from its groups to get a full set of detected genes.-------------
genes <- unique(c(rownames(spe_713_FFPE_filtered), rownames(spe_709_FFPE_filtered), rownames(spe_710_filtered), rownames(spe_708_FFPE_filtered),
                  rownames(spe_460_filtered), rownames(spe_462_filtered), rownames(spe_463_filtered),
                  rownames(spe_713_FFPE_CA_filtered), rownames(spe_709_FFPE_CA_filtered),
                  rownames(spe_709_FF_filtered), rownames(spe_708_FF_filtered), rownames(spe_545_filtered), rownames(spe_544_filtered),
                  rownames(spe_167_filtered),rownames(spe_168_filtered),rownames(spe_173_filtered),rownames(spe_174_filtered),
                  rownames(spe_709_FF_CA_filtered), rownames(spe_713_FF_CA_filtered)))

#------------4. Make a Venn diagram of group-level detected genes (from 2).-------------
## b. 2-way for FFPE and FFPE CA
venn_FFPE <- genes_grouped[1:2] %>%
  ggvenn(
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 6, text_size=7
  )

## c. 2-way for OCT and OCT CA
venn_OCT <- genes_grouped[3:4] %>%
  ggvenn(
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 6, text_size=7
  )

#-------------------Make figure--------------------
probes <- readxl::read_xlsx(here("analysis","Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.probe_metadata.xlsx")) # FFPE, OCT_CA, FFPE_CA
probes_included <- probes[probes$included==TRUE,]

# oct <- data.frame(gene_name=OCT_genes[OCT_genes %in% probes_included$gene_name])
probe_based <- unique(c(FFPE_genes[FFPE_genes %in% probes_included$gene_name],
                        FFPE_CA_genes[FFPE_CA_genes %in% probes_included$gene_name],
                        OCT_CA_genes[OCT_CA_genes %in% probes_included$gene_name]))

list_input <- list(OCT=OCT_genes, OCT_CA = OCT_CA_genes, FFPE=FFPE_genes, FFPE_CA=FFPE_CA_genes)

library(ComplexUpset)

pdf(here("analysis","output","detected_genes_upset_10pct.pdf"), height=6, width=10)
upset(UpSetR::fromList(list_input), intersect=c("OCT","OCT_CA","FFPE","FFPE_CA"),
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
                                    'overall_sizes'=theme(axis.text.x=element_text(size=17), axis.title.x = element_text(size=17)))))
dev.off()

library(ggvenn)
x <- list(polyA_based = OCT_genes, probe_based = probe_based, probe_set = probes_included$gene_name)
# pdf(here("analysis","output","detected_genes_venn.pdf"), height=6, width=8)
ggvenn(x,
       fill_color = c("#0073C2FF", "#EFC000FF","#CD534CFF"),
       stroke_size = 0.5, set_name_size = 5, text_size=4)
# dev.off()

# pdf(here("analysis","output","detected_genes_venn_OCT_10pct.pdf"))
venn_OCT
# dev.off()

#-----------compile figure--------------
library(cowplot)
top_row <- plot_grid(venn_all, venn_FFPE, labels=c("a","c"), label_size=30, rel_widths = c(1.3,0.8))
bottom_row <- plot_grid(upset, venn_OCT, labels=c("b","d"), label_size=30, rel_widths = c(1.3,0.8))

# pdf(here("analysis","output","10pct_gene_overlap.pdf"), height=12, width=16)#, units="in",res=1000)
plot_grid(top_row, bottom_row, nrow=2)
# dev.off()


#--------------------------------------
# names of OCT-only genes
elements <- unique(unlist(list_input))
a <- UpSetR::fromList(list_input)
oct_only <- elements[as.numeric(rownames(a[(a$OCT==1 & a$OCT_CA==0 & a$FFPE==0 & a$FFPE_CA==0),]))]
is_mito <- str_detect(oct_only, "Mt")
table(is_mito)
is_ribo <- str_detect(oct_only, "Rp")
table(is_ribo)
is_mito_ribo <- str_detect(oct_only, "Mrp")
table(is_mito_ribo)

is_either <- c(oct_only[str_detect(oct_only, "Rp")], oct_only[str_detect(oct_only, "Mrp")])

oct_only[!(oct_only %in% is_either)] %in% probes_included$gene_name
table(oct_only[!(oct_only %in% is_either)] %in% probes_included$gene_name)
