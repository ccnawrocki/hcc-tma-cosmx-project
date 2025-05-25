##### Cell Type Validating #####
## Cole Nawrocki ##

# Environment
rm(list = ls())
.libPaths()

# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project-final/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"

### Packages ###
## Data Science
library(plyr) # Tools
library(magrittr) # Pipes
library(tidyverse) # Tools
library(data.table) # Fast reading
library(Matrix) # Sparse data
## Plotting
library(ComplexHeatmap) # Heatmaps
## Bioinformatics
library(Seurat) # Processing
library(biomaRt) # Genome queries
library(InSituType) # Cell-typing

# CosMx directory
cosmx_dir <- "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx"

# Custom Functions
source(paste(cosmx_dir, "analysis/functions/Spatial_Functions.R", sep = "/"))

# Data
tmas <- LoadSeuratRds("hcc_tmas_final.RDS")
celltype_cols <- c("plasmablast"="magenta",
                   "tumor"="green4",         
                   "macrophage"="yellow",
                   "endothelial"="red3", 
                   "fibroblast"="grey", 
                   "mast"="gold4",
                   "B-cell"="purple",
                   "neutrophil"="green", 
                   "cholangiocyte"="dodgerblue", 
                   "Treg"="pink",
                   "mDC"="yellow3",
                   "T CD8 naive"="orange",
                   "T CD4 memory"="brown2",
                   "monocyte"="darkblue",     
                   "T CD8 memory"="beige",
                   "NK"="tan",
                   "pDC"="orange3", 
                   "T CD4 naive"="turquoise", 
                   "erythrocyte"="pink4")  

pdf("hcc_tmas_celltyping_flowcell.pdf", width = 12, height = 8)
FlowCellCellPlot(.METADATA = tmas@meta.data, 
                 .cellcolor = "celltype", 
                 .fovlabelcolor = "black",
                 size = 0.5, 
                 .fovlabel = F) + 
  scale_fill_manual(values = celltype_cols) +
  theme_void() + 
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  facet_grid(.~slide, scales = "free")
dev.off()

pdf("hcc_tmas_celltyping_UMAP.pdf", width = 8, height = 8)
DimPlot(tmas, group.by = "celltype", raster = F) +
  scale_color_manual(values = celltype_cols) +
  theme_void() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()

# Printing the tma1 images to check the cell-typing
images_setup(imdir_path = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc_tma1/small_rgb_images")

pdf("hcc_tma1_IF_celltype_plots.pdf", width = 12, height = 8)
for (f in (unique(tmas@meta.data |> filter(slide == "hcc_tma1") %$% fov |> as.numeric()))) {
  ip <- image_plot(.fov = f) + ggtitle(paste("FOV", f, sep = " "))
  
  cp <- FlowCellCellPlot(.METADATA = tmas@meta.data |> filter(slide == "hcc_tma1"), 
                         .cellcolor = "celltype",
                         .fovs = f, 
                         .fovlabelcolor = "black",
                         size = 2, 
                         .fovlabel = F) + 
    scale_fill_manual(values = celltype_cols) +
    theme_void() + 
    theme(panel.background = element_rect(fill = "black")) +
    guides(fill = guide_legend(override.aes = list(size = 5))) + 
    coord_fixed(expand = F)
  
  print(ip|cp)
}
dev.off()

# Erythrocytes
# Checking with an image... looks good
ip <- image_plot(.fov = 78)

cp <- FlowCellCellPlot(.METADATA = tmas@meta.data |> filter(slide == "hcc_tma1"), 
                       .cellcolor = "celltype",
                       .fovs = 78, 
                       .fovlabelcolor = "black",
                       size = 2, 
                       .fovlabel = F) + 
  scale_fill_manual(values = celltype_cols) +
  theme_void() + 
  theme(panel.background = element_rect(fill = "black")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  coord_fixed(expand = F)

print(ip|cp)

# Printing the tma2 images to check the cell-typing
images_setup(imdir_path = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc_tma2/small_rgb_images")
pdf("hcc_tma2_IF_celltype_plots.pdf", width = 12, height = 8)
for (f in (unique(tmas@meta.data |> filter(slide == "hcc_tma2") %$% fov |> as.numeric()))) {
  ip <- image_plot(.fov = f) + ggtitle(paste("FOV", f, sep = " "))
  
  cp <- FlowCellCellPlot(.METADATA = tmas@meta.data |> filter(slide == "hcc_tma2"), 
                         .cellcolor = "celltype",
                         .fovs = f, 
                         .fovlabelcolor = "black",
                         size = 2, 
                         .fovlabel = F) + 
    scale_fill_manual(values = celltype_cols) +
    theme_void() + 
    theme(panel.background = element_rect(fill = "black")) +
    guides(fill = guide_legend(override.aes = list(size = 5))) + 
    coord_fixed(expand = F)
  
  print(ip|cp)
}
dev.off()

