# Spot swapping with SpotClean
# Author: Mei Du

library(qs)
# all_samples <- qs::qload("all_samples.qs")

library(SpotClean) # need Seurat v5

# Load raw data
sample <- c(460,462,463,544,545,708,709,173,174,167,168,708,709,710,713,709,713,709,713)
run <- c(rep("V11L12-124",3),
         rep("V11D13-369",4),
         rep("V12F28-091",4),
         rep("V12M14-014",4),
         rep("V42L05-065",2),
         rep("V42L05-390",2))
protocol <- c(rep("FFPE manual",3),
              rep("OCT manual",4),
              rep("OCT manual",4),
              rep("FFPE manual",4),
              rep("OCT CA",2),
              rep("FFPE CA",2))
experiment <- c(rep("FFPE manual (earlier)", 3),
                rep("OCT manual (WT)",4),
                rep("OCT manual (KOvsCTRL)",4),
                rep("FFPE manual (later)", 4),
                rep("CytAssist",4))

metadata <- as.data.frame(cbind(sample, run, protocol, experiment))

clean_spots <- function(x) {
  spe <- read10xRaw(paste0(files[str_detect(files,paste0(metadata[x,"run"],"/.*",metadata[x,"sample"]))],"/outs/raw_feature_bc_matrix"))

  slide_info <- read10xSlide(tissue_csv_file = paste0(files[str_detect(files,paste0(metadata[x,"run"],"/.*",metadata[x,"sample"]))], "/outs/spatial/tissue_positions.csv"),
                             tissue_img_file=paste0(files[str_detect(files,paste0(metadata[x,"run"],"/.*",metadata[x,"sample"]))], "/outs/spatial/tissue_lowres_image.png"),
                             scale_factor_file=paste0(files[str_detect(files,paste0(metadata[x,"run"],"/.*",metadata[x,"sample"]))], "/outs/spatial/scalefactors_json.json"))
  slide_obj <- createSlide(spe, slide_info)
  decont_obj <- spotclean(slide_obj, maxit=10, candidate_radius = 20)
}

sc_list <- list()
for (i in 1:nrow(metadata)) {
  sc_list[[i]] <- clean_spots(i)
}
saveRDS(sc_list, "output/RDS/spotclean.RDS")

# Output metrics

# per-spot contamination rate
metadata(decont_obj)$contamination_rate

bleed <- data.frame(bleeding_rate=0)
cont <- data.frame(mean_spot_contamination=0)
for (i in 1:length(sc_list)) {
  bleed[i,"bleeding_rate"] <- metadata(sc_list[[i]])$bleeding_rate
  cont[i, "mean_spot_contamination"] <- mean(metadata(sc_list[[i]])$contamination_rate)
}

cbind(metadata,bleed,cont)

