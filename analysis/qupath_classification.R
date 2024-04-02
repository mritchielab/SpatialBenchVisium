# Visualising QuPath classification
# Author: Mei Du

library(here)
library(dplyr)
library(ggplot2)
library(SpatialExperiment)
library(patchwork)

##--------------------------qupath---------------------------##
# install.packages("sf",configure.args= "--with-gdal-config=/stornext/System/data/tools/gdal/gdal-2.4.4/bin/gdal-config")
library(sf)

file_709 <- st_read("FFPE_709.geojson") %>% st_transform(., src=4236 )
file_713 <- st_read("FFPE_713.geojson") %>% st_transform(., src=4236 )
file_709$classification <- c("Red Pulp","White Pulp")
file_713$classification <- c("Red Pulp","White Pulp")


plot_709 <- ggplot() +
  geom_sf(data = file_709, aes(fill=classification)) +
  scale_fill_manual(values = c("Red Pulp" = "red", "White Pulp" = "white")) +
  theme_minimal() +
  theme(legend.position="none")
plot_713 <- ggplot() +
  geom_sf(data = file_713, aes(fill=classification)) +
  scale_fill_manual(values = c("Red Pulp" = "red", "White Pulp" = "white")) +
  theme_minimal()

# manually increase width
# pdf(here("analysis","output","qupath_classification.pdf"))
plot_709 + plot_713
# dev.off()
