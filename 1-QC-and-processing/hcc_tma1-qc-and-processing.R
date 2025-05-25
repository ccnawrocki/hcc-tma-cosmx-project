rm(list = ls())

### Packages ###
## Data Science
library(magrittr) # Pipes
library(tidyverse) # Tools
library(data.table) # Fast reading
library(Matrix) # Sparse data
library(arrow) # Lazy loading
## Spatial Data
library(sf) # Segmentations and shapes
## Plotting
library(patchwork) # Combining plots
library(ggridges) # Ridge plots
library(ggsankey) # Sankey plots
library(circlize) # Color ramps
library(ComplexHeatmap) # Heatmaps
library(ggpubr) # Image plots
## Bioinformatics
library(Seurat) # Processing

# CosMx directory
cosmx_dir <- "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx"

# Custom Functions
source(paste(cosmx_dir, "analysis/functions/Spatial_Functions.R", sep = "/"))

# Environment Path
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project-final/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"   

## Data Overview ---------------------------------------------------------------
# Examining the flatFiles
data_dir <- paste(cosmx_dir, "franses-hcc-data/hcc_tma1", sep = "/") # data directory
list.files(paste(data_dir, "flatFiles/hcc_tma1", sep = "/")) # There should be 5 files here

# Setting up the tx file for easy queries (this is a custom function)
# Creating an arrow dataset to query only the tx needed at a given time, which saves memory
#SetupCosMxTx(.DIR = paste(data_dir, "flatFiles", sep = "/"), .NAME = "hcc_tma1") 

# Reading the data into memory (this is a custom function)
tma1 <- ReadCosMxFF(.DIR = paste(data_dir, "flatFiles", sep = "/"), .NAME = "hcc_tma1")

# Setting up the segmentation for easier plotting (this is a custom function)
tma1 <- SegmentationToMetadata(.SLIDE = tma1)

# Making the segmentation readable for QuPath (this is a custom function)
#MakeGeoJsons(.SLIDE = tma1, .out_dir = paste(data_dir, "geojson_segmentations", sep = "/"))

# Adding metadata
coremap <- readxl::read_excel("1-QC-and-processing/hcc_tma_maps.xlsx", 
                              sheet = "tma1-pathology")
tma1$metadata$fov_type <- plyr::mapvalues(x = tma1$metadata$fov, from = coremap$FOV, to = coremap$Inference)
tma1$metadata$fov_type %<>% as.factor() # Making the new column a factor
tma1$metadata$core <- plyr::mapvalues(x = tma1$metadata$fov, from = coremap$FOV, to = coremap$Comments)
fov_type_map <- tma1$metadata |> dplyr::group_by(fov) |> summarise(fov_type = unique(fov_type))
tma1$metadata$core %<>% as.factor() # Making the new column a factor
pt_map <- openxlsx::read.xlsx("1-QC-and-processing/hcc_tma_maps.xlsx", sheet = "tma1-map", cols = 6:8)
tma1$metadata$patient <- plyr::mapvalues(x = tma1$metadata$core, from = pt_map$core, to = pt_map$patient) |> as.factor()
# The following `from` values were not present in `x`: c9, c11, c13, c23, c24

# Color map for plotting the fov types
fov_type_cmap <- c("T"="green4", "EFT"="lightgreen", "EFN"="cyan", "N"="dodgerblue", "S"="gold", "BNT"="red", "EFBNT"="salmon", "Kidney"="purple")

# Plots of core and patient... looks good (custom function)
FlowCellCellPlot(tma1$metadata, .cellcolor = "core", size = 1, .fovlabel = F) + 
  ggprism::scale_fill_prism() + 
  theme_void() + 
  coord_fixed() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
FlowCellCellPlot(tma1$metadata, .cellcolor = "patient", size = 1, .fovlabel = F) + 
  ggprism::scale_fill_prism() + 
  theme_void() + 
  coord_fixed() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))

# Plotting depth-related statistics on the flowcell, using a custom function
sop1 <- FlowCellCellPlot(.METADATA = tma1$metadata, .cellcolor = "nCount_RNA", .fovlabel = F)  + 
  scale_fill_gradient2(low="beige", mid="beige", high="dodgerblue", 
                       limits = quantile(tma1$metadata$nCount_RNA, probs = c(0, 0.975)), 
                       oob = scales::squish) + 
  theme_void() + 
  coord_fixed() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
sop1 <- FlowCellFOVPlot(.METADATA = tma1$metadata, .fovcolor = "fov_type", .fovlabel = F, .cellplot = sop1) + 
  scale_color_manual(values = fov_type_cmap) + coord_fixed()
sop1

sop2 <- FlowCellCellPlot(.METADATA = tma1$metadata, .cellcolor = "nFeature_RNA", .fovlabel = F) + 
  scale_fill_gradient2(low="beige", mid="beige", high="red", 
                       limits = quantile(tma1$metadata$nFeature_RNA, probs = c(0, 0.975)), 
                       oob = scales::squish) + 
  theme_void() + 
  coord_fixed() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
sop2 <- FlowCellFOVPlot(.METADATA = tma1$metadata, .fovcolor = "fov_type", .fovlabel = F, .cellplot = sop2) + 
  scale_color_manual(values = fov_type_cmap) + coord_fixed()
sop2

# Plotting background-related statistics on the flowcell, using a custom function
sop3 <- FlowCellCellPlot(.METADATA = tma1$metadata, .cellcolor = "nCount_negprobes", .fovlabel = F)  + 
  scale_fill_gradient2(low="beige", mid="beige", high="magenta", 
                       limits = quantile(tma1$metadata$nCount_negprobes, probs = c(0, 0.975)), 
                       oob = scales::squish) + 
  theme_void() + 
  coord_fixed() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
sop3 <- FlowCellFOVPlot(.METADATA = tma1$metadata, .fovcolor = "fov_type", .fovlabel = F, .cellplot = sop3) + 
  scale_color_manual(values = fov_type_cmap) + coord_fixed()
sop3

sop4 <- FlowCellCellPlot(.METADATA = tma1$metadata, .cellcolor = "nCount_falsecode", .fovlabel = F) + 
  scale_fill_gradient2(low="beige", mid="beige", high="orange", 
                       limits = quantile(tma1$metadata$nCount_falsecode, probs = c(0, 0.975)), 
                       oob = scales::squish) + 
  theme_void() + 
  coord_fixed() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
sop4 <- FlowCellFOVPlot(.METADATA = tma1$metadata, .fovcolor = "fov_type", .fovlabel = F, .cellplot = sop4) + 
  scale_color_manual(values = fov_type_cmap) + coord_fixed()
sop4

# Calculating total counts and genes
tma1$metadata$total_counts <- tma1$counts |> colSums()
tma1$metadata$total_features <- (tma1$counts > 0) |> colSums()

## Cell-Level QC ---------------------------------------------------------------
# We will filter cells on total RNA counts, cell area, proportion of counts from background, and log10(counts/positive gene).
before <- ncol(tma1$counts)

# Identifying negative probes and system controls
idx <- grepl(rownames(tma1$counts), pattern = "Negative")
idx2 <- grepl(rownames(tma1$counts), pattern = "SystemControl")
negs <- tma1$counts[idx,]
scontrs <- tma1$counts[idx2,]
rnacts <- tma1$counts[!(idx|idx2),]

# Plotting the distributions of these variables first
tma1$metadata$fov %<>% as.factor() # Converting the fov column to factor for plotting purposes

# total counts
cqcp1 <- ggplot() +
  geom_density_ridges_gradient(data = tma1$metadata, mapping = aes(x = nCount_RNA, y = fov, fill = fov_type)) +  
  facet_wrap(.~fov_type, scales = "free") +
  geom_vline(xintercept = c(20, 2000), color = "red2") + 
  theme_classic() + scale_fill_manual(values = fov_type_cmap)

# Area
cqcp2 <- ggplot() +
  geom_density_ridges_gradient(data = tma1$metadata, mapping = aes(x = Area.um2, y = fov, fill = fov_type)) +  
  facet_wrap(.~fov_type, scales = "free") +
  geom_vline(xintercept = c(25, 1000), color = "cyan2") + 
  theme_classic() + scale_fill_manual(values = fov_type_cmap)

# Proportion of counts from background
tma1$metadata$prop_counts_negative <- tma1$metadata$nCount_negprobes/(tma1$metadata$nCount_RNA + tma1$metadata$nCount_negprobes)
cqcp3 <- ggplot() +
  geom_density_ridges_gradient(data = tma1$metadata, mapping = aes(x = prop_counts_negative, y = fov, fill = fov_type)) +  
  facet_wrap(.~fov_type, scales = "free") +
  geom_vline(xintercept = c(0.1), color = "magenta2") + 
  theme_classic() + scale_fill_manual(values = fov_type_cmap)

# Complexity score: log10(counts per positive gene)
tma1$metadata$log10counts_per_gene <- (tma1$metadata$nCount_RNA/tma1$metadata$nFeature_RNA) |> log10()
cqcp4 <- ggplot() +
  geom_density_ridges_gradient(data = tma1$metadata, mapping = aes(x = log10counts_per_gene, y = fov, fill = fov_type)) +  
  facet_wrap(.~fov_type, scales = "free") +
  geom_vline(xintercept = 2, color = "gold2") + 
  theme_classic() + scale_fill_manual(values = fov_type_cmap)

# Plotting (uses patchwork)
cqcp1/cqcp2/cqcp3/cqcp4

# Plotting the proportion of each FOV that will be removed
# This function will return if each cell passes on each of the four filtering metrics 
cellpassqc <- function(.metadata) {
  pass_counts_qc <- ((.metadata$nCount_RNA > 20) & (.metadata$nCount_RNA < 2000))
  pass_area_qc <- ((.metadata$Area.um2 > 25) & (.metadata$Area.um2 < 1000))
  pass_bg_qc <- (.metadata$prop_counts_negative < 0.1)
  pass_log10cpg_qc <- (.metadata$log10counts_per_gene < 2)
  return(list("pass_counts_qc" = pass_counts_qc, 
              "pass_area_qc" = pass_area_qc, 
              "pass_bg_qc" = pass_bg_qc, 
              "pass_log10cpg_qc" = pass_log10cpg_qc))
}

# Adding the results to the metadata
qcresults <- cellpassqc(.metadata = tma1$metadata)
tma1$metadata %<>% cbind(qcresults)
tma1$metadata %<>% mutate(passes_qc = ((pass_counts_qc+pass_area_qc+pass_bg_qc+pass_log10cpg_qc) == 4)) 

# Transforming the data
qcd <- tma1$metadata |> 
  pivot_longer(cols = c(pass_counts_qc, pass_area_qc, pass_bg_qc, pass_log10cpg_qc, passes_qc), names_to = "metric", values_to = "pass") |> 
  dplyr::group_by(fov, metric, pass) |> 
  tally() |> 
  pivot_wider(names_from = pass, values_from = n, values_fill = 0) |> 
  mutate(fail_rate = `FALSE`/(`FALSE`+`TRUE`))
qcd$fov_type <- plyr::mapvalues(x = qcd$fov, from = fov_type_map$fov, to = as.character(fov_type_map$fov_type))

# Plotting
cqcp5 <- ggplot(data = qcd) + 
  geom_point(mapping = aes(x = fov, y = fail_rate, color = fov_type, label = fov)) +
  ggrepel::geom_text_repel(mapping = aes(x = fov, y = fail_rate, color = fov_type, label = fov)) + 
  facet_wrap(.~metric, nrow = 3, ncol = 2) + 
  scale_color_manual(values = fov_type_cmap) +
  theme_classic() + 
  theme(axis.text.x = element_blank())

# Up-close view of the segmentation of the fovs with fail rate > 25%, which will be removed
toremove <- qcd |> filter(metric == "passes_qc" & fail_rate > 0.25) %$%fov
cqcp6 <- SegmentationPlot(.METADATA = tma1$metadata, .fovs = toremove, .cellcolor = "passes_qc", .outlinecolor = "black", ncol = 2) + 
  scale_fill_manual(values = c("black", "grey")) + 
  theme_classic()

cqcp5
cqcp6

# Actually doing the filtering
tma1$metadata %<>% filter(passes_qc == T)
tma1$metadata %<>% filter(!(fov %in% toremove))
tma1$counts <- tma1$counts[,rownames(tma1$metadata)]

after <- ncol(tma1$counts)
after/before
# [1] 0.9158034

## Gene-level QC ---------------------------------------------------------------
# A gene is flagged if its total expression is not greater than the mean negative probe 
# total expression plus two standard deviations of the total negative probe expressions. 
# Simply put, for gene X, if (total counts X) < (mean total counts negative probes) + 
# 2*(standard deviation total counts negative probes), then gene X will be flagged. This
# gene should not be used as a marker gene based on positivity alone. 

# Getting negative probes again
negs <- tma1$counts[idx,]

# Getting the threshold
mean_neg <- negs |> rowSums() |> mean()
sd_neg <- negs |> rowSums() |> sd()
thresh <- mean_neg + 2*sd_neg

# Deciding if each gene is above background
above_bg <- (tma1$counts[!(idx|idx2),] |> rowSums() > thresh)
above_bg |> sum()
# [1] 589

# Identifying genes whose proportion of counts outside of cells are not significantly less than 0.5 as a precautionary measure
cellcompcts <- tma1$tx |> dplyr::group_by(fov, target, CellComp) |> tally() |> collect() |> pivot_wider(names_from = CellComp, values_from = n, values_fill = 0)
cellcompcts <- cellcompcts |> mutate(incells = Membrane+Cytoplasm+Nuclear)
cellcompcts <- cellcompcts |> mutate(background = None)
cellcompcts$bg_prop <- cellcompcts$background/(cellcompcts$incells + cellcompcts$background)
cellcompcts <- cellcompcts |> dplyr::group_by(target) |> dplyr::summarise(meanprop = mean(bg_prop), sdprop = sd(bg_prop))
cellcompcts$lower <- cellcompcts$meanprop-(qt(0.975, df = n_distinct(tma1$metadata$fov))-1)*(cellcompcts$sdprop/sqrt(n_distinct(tma1$metadata$fov)))
cellcompcts$upper <- cellcompcts$meanprop+(qt(0.975, df = n_distinct(tma1$metadata$fov))-1)*(cellcompcts$sdprop/sqrt(n_distinct(tma1$metadata$fov)))
idx3 <- cellcompcts[cellcompcts$upper > 0.5, ]$target
outgenes <- idx3
idx3 <- rownames(tma1$counts) %in% idx3

# Subsetting the counts, based on the work from above
tma1$counts <- tma1$counts[!(idx|idx2|idx3),]

# Making sure no genes have 0 counts
zerocountgenes <- (tma1$counts |> rowSums() == 0) |> sum()
outgenes
# character(0)
zerocountgenes
# [1] 0

## Processing with Seurat ------------------------------------------------------
# Using the Seurat v3 object, since it is easiest to use
# Standard Seurat workflow
# - Normalization: Log2 counts per 1000 total counts per cell normalization
# - Variable genes: Use all genes as variable genes, identified via dispersion
# - Scaling: Only variable genes were scaled, without regressing any variables out
# - Dimensional reduction: PCA was run, using the variable genes, and then the top 25 PCs were used for UMAP formation
options(Seurat.object.assay.version = "v3")
library(parallel)
options(mc.cores = 8)

# Creating the object
# *Cannot have the sf polygon column in the metadata when the object is created*
# Best to do for control separately
metadata <- dplyr::select(tma1$metadata, c(1:17, 26, 27, 57, 58, 60:71))
tma1so <- CreateSeuratObject(counts = tma1$counts, meta.data = metadata, project = "hcc_tma1")
control <- subset(tma1so, subset = fov_type == "Kidney")
tma1so <- subset(tma1so, subset = fov_type != "Kidney")

# Normalizing
scale_factor <- 1000
tma1so <- NormalizeData(tma1so, normalization.method = "LogNormalize", scale.factor = scale_factor)
control <- NormalizeData(control, normalization.method = "LogNormalize", scale.factor = scale_factor)

# Variable Features
tma1so <- FindVariableFeatures(tma1so, selection.method = "disp", nfeatures = nrow(tma1so))
control <- FindVariableFeatures(control, selection.method = "disp", nfeatures = nrow(control))

# Plotting the variable features
tma1so@assays$RNA@meta.features$gene <- rownames(tma1so@assays$RNA@meta.features) # Adding the gene names to the gene metadata
tma1so@assays$RNA@meta.features$below_bg <- !above_bg[rownames(tma1so)] # Adding whether or not a gene is above background, as defined above
tma1so@assays$RNA@meta.features$mvp.variable %<>% as.logical()
tma1so@assays$RNA@meta.features$below_bg %<>% as.logical()
control@assays$RNA@meta.features$gene <- rownames(control@assays$RNA@meta.features) # Adding the gene names to the gene metadata
control@assays$RNA@meta.features$below_bg <- !above_bg[rownames(control)] # Adding whether or not a gene is above background, as defined above
control@assays$RNA@meta.features$mvp.variable %<>% as.logical()
control@assays$RNA@meta.features$below_bg %<>% as.logical()

# Making the plot
vfp <- ggplot(data = tma1so@assays$RNA@meta.features) + 
  geom_point(mapping = aes(x=mvp.mean, y=mvp.dispersion, color = mvp.variable, shape = below_bg)) + 
  scale_color_manual(values = c("black", "red")) + 
  ggrepel::geom_text_repel(data = tma1so@assays$RNA@meta.features |> filter(mvp.variable == T), 
                           mapping = aes(x = mvp.mean, y = mvp.dispersion, label = gene))
vfp
vfp <- ggplot(data = control@assays$RNA@meta.features) + 
  geom_point(mapping = aes(x=mvp.mean, y=mvp.dispersion, color = mvp.variable, shape = below_bg)) + 
  scale_color_manual(values = c("black", "red")) + 
  ggrepel::geom_text_repel(data = control@assays$RNA@meta.features |> filter(mvp.variable == T), 
                           mapping = aes(x = mvp.mean, y = mvp.dispersion, label = gene))
vfp

# Scaling
tma1so <- ScaleData(tma1so, features = VariableFeatures(tma1so))
control <- ScaleData(control, features = VariableFeatures(control))

# PCA
# Using 25 PCs will be plenty
tma1so <- RunPCA(tma1so, features = VariableFeatures(tma1so), npcs = 50)
ElbowPlot(tma1so, ndims = 50)
# Using 15 PCs will be plenty
control <- RunPCA(control, features = VariableFeatures(control), npcs = 50)
ElbowPlot(control, ndims = 50)

# UMAP
tma1so <- RunUMAP(tma1so, dims = 1:25)
DimPlot(tma1so, group.by = "fov_type", raster = F) + scale_color_manual(values = fov_type_cmap) + theme_void()
DimPlot(tma1so, group.by = "core", raster = F) + ggprism::scale_color_prism() + theme_void()
DimPlot(tma1so, group.by = "patient", raster = F, split.by = "patient") + ggprism::scale_color_prism() + theme_void()
FeaturePlot(tma1so, features = "nCount_RNA") + scale_color_viridis_c(option = "turbo") + theme_void()
control <- RunUMAP(control, dims = 1:15)
DimPlot(control, group.by = "fov_type", raster = F) + scale_color_manual(values = fov_type_cmap) + theme_void()
DimPlot(control, group.by = "core", raster = F) + ggprism::scale_color_prism() + theme_void()
DimPlot(control, group.by = "patient", raster = F, split.by = "patient") + ggprism::scale_color_prism() + theme_void()
FeaturePlot(control, features = "nCount_RNA") + scale_color_viridis_c(option = "turbo") + theme_void()

tma1so[["NP"]] <- CreateAssayObject(counts = negs[,colnames(tma1so)])
tma1so[["FC"]] <- CreateAssayObject(counts = scontrs[,colnames(tma1so)])
SaveSeuratRds(object = tma1so, file = "hcc_tma1_final.RDS")

control[["NP"]] <- CreateAssayObject(counts = negs[,colnames(control)])
control[["FC"]] <- CreateAssayObject(counts = scontrs[,colnames(control)])
SaveSeuratRds(object = control, file = "hcc_tma1_control_final.RDS")

## Session ---------------------------------------------------------------------
# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS 15.1.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] parallel  grid      stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] InSituType_2.0        biomaRt_2.62.0        Seurat_5.1.0          SeuratObject_5.0.2    sp_2.1-4              ggpubr_0.6.0         
# [7] ComplexHeatmap_2.22.0 circlize_0.4.16       ggsankey_0.0.99999    ggridges_0.5.6        patchwork_1.3.0       sf_1.0-19            
# [13] arrow_18.1.0.1        Matrix_1.7-1          data.table_1.16.4     lubridate_1.9.4       forcats_1.0.0         stringr_1.5.1        
# [19] dplyr_1.1.4           purrr_1.0.2           readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1        
# [25] tidyverse_2.0.0       magrittr_2.0.3       
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.5                    matrixStats_1.4.1           spatstat.sparse_3.1-0       httr_1.4.7                 
# [5] RColorBrewer_1.1-3          doParallel_1.0.17           tools_4.4.1                 sctransform_0.4.1          
# [9] backports_1.5.0             R6_2.5.1                    lazyeval_0.2.2              uwot_0.2.2                 
# [13] GetoptLong_1.0.5            withr_3.0.2                 prettyunits_1.2.0           gridExtra_2.3              
# [17] rematch_2.0.0               progressr_0.15.0            cli_3.6.3                   Biobase_2.66.0             
# [21] spatstat.explore_3.3-3      fastDummies_1.7.4           labeling_0.4.3              spatstat.data_3.1-4        
# [25] proxy_0.4-27                pbapply_1.7-2               askpass_1.2.1               R.utils_2.12.3             
# [29] parallelly_1.39.0           readxl_1.4.3                rstudioapi_0.17.1           RSQLite_2.3.9              
# [33] generics_0.1.3              shape_1.4.6.1               ica_1.0-3                   spatstat.random_3.3-2      
# [37] zip_2.3.1                   car_3.1-3                   S4Vectors_0.44.0            abind_1.4-8                
# [41] R.methodsS3_1.8.2           lifecycle_1.0.4             carData_3.0-5               SummarizedExperiment_1.36.0
# [45] SparseArray_1.6.0           BiocFileCache_2.14.0        Rtsne_0.17                  blob_1.2.4                 
# [49] promises_1.3.0              crayon_1.5.3                miniUI_0.1.1.1              lattice_0.22-6             
# [53] cowplot_1.1.3               KEGGREST_1.46.0             pillar_1.10.1               GenomicRanges_1.58.0       
# [57] rjson_0.2.23                future.apply_1.11.3         codetools_0.2-20            leiden_0.4.3.1             
# [61] glue_1.8.0                  spatstat.univar_3.1-1       vctrs_0.6.5                 png_0.1-8                  
# [65] spam_2.11-0                 cellranger_1.1.0            gtable_0.3.6                assertthat_0.2.1           
# [69] cachem_1.1.0                openxlsx_4.2.7.1            S4Arrays_1.6.0              mime_0.12                  
# [73] survival_3.7-0              SingleCellExperiment_1.28.1 iterators_1.0.14            units_0.8-5                
# [77] fitdistrplus_1.2-1          ROCR_1.0-11                 lsa_0.73.3                  nlme_3.1-166               
# [81] bit64_4.5.2                 progress_1.2.3              filelock_1.0.3              RcppAnnoy_0.0.22           
# [85] GenomeInfoDb_1.42.1         SnowballC_0.7.1             irlba_2.3.5.1               KernSmooth_2.23-24         
# [89] colorspace_2.1-1            BiocGenerics_0.52.0         DBI_1.2.3                   tidyselect_1.2.1           
# [93] bit_4.5.0                   compiler_4.4.1              curl_6.1.0                  httr2_1.0.7                
# [97] xml2_1.3.6                  DelayedArray_0.32.0         plotly_4.10.4               scales_1.3.0               
# [101] classInt_0.4-10             lmtest_0.9-40               rappdirs_0.3.3              digest_0.6.37              
# [105] goftest_1.2-3               spatstat.utils_3.1-2        XVector_0.46.0              htmltools_0.5.8.1          
# [109] pkgconfig_2.0.3             umap_0.2.10.0               MatrixGenerics_1.18.0       dbplyr_2.5.0               
# [113] fastmap_1.2.0               rlang_1.1.4                 GlobalOptions_0.1.2         htmlwidgets_1.6.4          
# [117] UCSC.utils_1.2.0            shiny_1.9.1                 farver_2.1.2                zoo_1.8-12                 
# [121] jsonlite_1.8.9              mclust_6.1.1                R.oo_1.27.0                 Formula_1.2-5              
# [125] GenomeInfoDbData_1.2.13     dotCall64_1.2               munsell_0.5.1               Rcpp_1.0.13-1              
# [129] reticulate_1.40.0           stringi_1.8.4               zlibbioc_1.52.0             MASS_7.3-61                
# [133] plyr_1.8.9                  listenv_0.9.1               ggrepel_0.9.6               deldir_2.0-4               
# [137] Biostrings_2.74.1           splines_4.4.1               tensor_1.5                  hms_1.1.3                  
# [141] igraph_2.1.1                spatstat.geom_3.3-4         ggsignif_0.6.4              RcppHNSW_0.6.0             
# [145] reshape2_1.4.4              stats4_4.4.1                renv_1.0.11                 BiocManager_1.30.25        
# [149] ggprism_1.0.5               tzdb_0.4.0                  foreach_1.5.2               httpuv_1.6.15              
# [153] RANN_2.6.2                  openssl_2.3.1               polyclip_1.10-7             future_1.34.0              
# [157] clue_0.3-66                 scattermore_1.2             broom_1.0.7                 xtable_1.8-4               
# [161] e1071_1.7-16                RSpectra_0.16-2             rstatix_0.7.2               later_1.3.2                
# [165] viridisLite_0.4.2           class_7.3-22                memoise_2.0.1               AnnotationDbi_1.68.0       
# [169] IRanges_2.40.1              cluster_2.1.6               timechange_0.3.0            globals_0.16.3             

