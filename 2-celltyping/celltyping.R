##### Cell-typing #####
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

# Data 
tma1so <- LoadSeuratRds("hcc_tma1_final.RDS")
tma2so <- LoadSeuratRds("hcc_tma2_final.RDS")

# TMA1 -------------------------------------------------------------------------
# Identifying tumor
tma1so <- FindNeighbors(object = tma1so, dims = 1:25)
tma1so <- FindClusters(tma1so, algorithm = 1, resolution = seq(0.5, 1.5, 0.1))
DimPlot(tma1so, raster = F, order = T,  label = T, group.by = "RNA_snn_res.0.5") + NoLegend()

FeaturePlot(object = tma1so, features = "Mean.PanCK", order = T, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma1so, features = "Mean.CD45", order = F, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma1so, features = "Mean.CD68_CK8_18", order = T, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma1so, features = "Area.um2", order = T, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma1so, features = "APOA1", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma1so, features = "GC", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma1so, features = "VTN", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma1so, features = "ARG1", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma1so, features = "FGG", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma1so, features = "GPC3", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()

tma1so <- AddModuleScore(tma1so, features = list(c("APOA1", "GC", "VTN", "ARG1", "FGG", "GPC3")), nbin = 5, name = "tumor.score")
FeaturePlot(tma1so, features = "tumor.score1", order = T) + scale_color_viridis_c()

suspected_tumor <- (tma1so$RNA_snn_res.0.5 %in% c("4", "0", "6", "9", "10", "15", "13", "12", "5"))
tma1so$suspected_tumor <- NA
tma1so$suspected_tumor[suspected_tumor] <- "tumor"
DimPlot(tma1so, group.by = "suspected_tumor", raster = F, order = T) 

# Identifying Erythrocytes
# Genes with three blue bits in their barcodes should all be expressed in them
probes <- read.table(file = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/MGH_Franses_panel.txt", 
                     sep = "\t", 
                     header = T)
probes$Barcode[probes$Barcode |> grepl(pattern = ".*BB.*BB.*BB.*") |> which()]
bluetargets <- probes$DisplayName[probes$Barcode |> grepl(pattern = ".*BB.*BB.*BB.*") |> which()] |> intersect(rownames(tma1so))
tma1so <- AddModuleScore(tma1so, features = list(bluetargets), name = "Autofluorescence.Score", nbin = 5)
FeaturePlot(tma1so, features = "Autofluorescence.Score1", raster = F, order = F) + scale_color_viridis_c()

suspected_rbc <- (tma1so$RNA_snn_res.0.5 == "11")
tma1so$suspected_rbc <- NA
tma1so$suspected_rbc[suspected_rbc] <- "rbc"
DimPlot(tma1so, group.by = "suspected_rbc", raster = F, order = T) 

# Loading in reference dataset
cosmx_dir <- "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx"
refcts <- fread(paste(cosmx_dir, "reference-datasets/human-protein-atlas/rna_single_cell_read_count/liver/read_count.tsv", sep = "/"))
refmeta <- fread(paste(cosmx_dir, "reference-datasets/human-protein-atlas/rna_single_cell_read_count/liver/cell_data.tsv", sep = "/"))
refclstinfo <- fread(paste(cosmx_dir, "reference-datasets/human-protein-atlas/rna_single_cell_cluster_description.tsv", sep = "/"))

# Cleaning it 
refcts <- refcts[-1,]
colnames(refcts) <- c("gene", paste("cell", 1:(ncol(refcts)-1), sep = "_"))
refcts <- as.data.frame(refcts)

ensembles <- refcts$gene
human <- useEnsembl(mirror = "useast", biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_key <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id", values = ensembles, mart=human)
refcts$gene <- mapvalues(x = refcts[[1]], from = gene_key$ensembl_gene_id, to = gene_key$hgnc_symbol)
refcts <- dplyr::filter(refcts, gene %in% rownames(tma1so))
refcts %<>% tibble::column_to_rownames("gene")

refmeta$cell_id %<>% paste("cell", ., sep = "_")
refclstinfo$Cluster %<>% gsub("c-", "", .) %<>% as.numeric()
refmeta$celltype <- mapvalues(x = refmeta$cluster, from = refclstinfo |> dplyr::filter(Tissue == "Liver") %$% Cluster, to = refclstinfo |> dplyr::filter(Tissue == "Liver") %$% `Cell type`)

# Getting reference profiles
mm <- model.matrix(~0+celltype, data = refmeta)
aggrefcts <- as.matrix(refcts) %*% mm
cellsums <- refmeta |> dplyr::group_by(celltype) |> dplyr::tally() %$% n
meanrefcts <- sweep(aggrefcts, MARGIN = 2, STATS = cellsums, FUN = "/")
colnames(meanrefcts) %<>% gsub("celltype", "", .)

refcts <- refcts |> t() |> as.data.frame()
sdrefcts <- refcts |> split(refmeta$celltype) |> purrr::map(as.matrix) |> purrr::map(matrixStats::colSds)
sdrefcts <- dplyr::bind_rows(sdrefcts, .id = "celltype") |> tibble::column_to_rownames("celltype") |> t()

# Defining cohort data
tma1so$sdimx <- tma1so$CenterX_global_px*(0.51/4256)
tma1so$sdimy <- tma1so$CenterY_global_px*(0.51/4256)

cidx <- (is.na(tma1so$suspected_tumor) & is.na(tma1so$suspected_rbc))
neighbors <- InSituCor:::nearestNeighborGraph(x = tma1so$sdimx, y = tma1so$sdimy, N = 50, subset = tma1so$core) 
neighborexpression <- InSituCor:::get_neighborhood_expression(counts = tma1so@assays$RNA@counts |> t(), neighbors = neighbors) 
neighborhoodPCs <- irlba::prcomp_irlba(neighborexpression, n = 10)$x

cohortdata <- tma1so@meta.data[cidx,] |> 
  dplyr::select(Area, Mean.PanCK, Mean.CD45, Mean.CD68_CK8_18, Mean.DAPI, Mean.CD298_B2M)
cohort <- InSituType::fastCohorting(mat = cbind(cohortdata, neighborhoodPCs[cidx,]))

neg <- tma1so@assays$NP@counts

# Running the function
sup <- insitutype(
  x = tma1so@assays$RNA@counts[,cidx] |> t(),
  neg = neg[,cidx] |> colMeans(),
  assay_type = "RNA", 
  n_clusts=0, 
  anchors = NULL, 
  reference_profiles = meanrefcts,
  reference_sds = sdrefcts,
  cohort = cohort,
  update_reference_profiles = T, 
  rescale = T, 
  refit = F
) 

tma1so$insitutype_sup_cluster <- NA
tma1so$insitutype_sup_cluster[cidx] <- sup$clust

tma1so$insitutype_sup_prob <- NA
tma1so$insitutype_sup_prob[cidx] <- sup$prob

DimPlot(tma1so, group.by = "insitutype_sup_cluster", raster = F, order = F)

data("ioprofiles")

cidx2 <- (!is.na(tma1so$insitutype_sup_cluster)) & (tma1so$insitutype_sup_cluster != "Cholangiocytes")

cohortdata2 <- tma1so@meta.data[cidx2,] |> 
  dplyr::select(Area, Mean.PanCK, Mean.CD45, Mean.CD68_CK8_18, Mean.DAPI, Mean.CD298_B2M)
cohort2 <- InSituType::fastCohorting(mat = cbind(cohortdata2, neighborhoodPCs[cidx2,]))

source("2-celltyping/getSubtypingGenes.R") # A script from NanoString
safegenes <- findSafeGenes(
  counts = t(tma1so@assays$RNA$counts),
  xy = tma1so@meta.data |> dplyr::select(sdimx, sdimy) |> as.matrix(),
  ismycelltype = (cidx2),
  tissue = tma1so$core,
  self_vs_neighbor_threshold = 1.75
)$safegenes           

sup2 <- insitutype(
  x = tma1so@assays$RNA@counts[safegenes,cidx2] |> t(),
  neg = neg[,cidx2] |> colMeans(),
  assay_type = "RNA", 
  n_clusts=0, 
  anchors = NULL, 
  reference_profiles = ioprofiles,
  cohort = cohort2,
  update_reference_profiles = T, 
  rescale = T, 
  refit = F
) 

tma1so$insitutype_sup_cluster2 <- NA
tma1so$insitutype_sup_cluster2[cidx2] <- sup2$clust

tma1so$insitutype_sup_prob2 <- NA
tma1so$insitutype_sup_prob2[cidx2] <- sup2$prob

DimPlot(tma1so, group.by = "insitutype_sup_cluster2", raster = F, order = F) + ggprism::scale_color_prism()

tma1so$celltype <- tma1so$insitutype_sup_cluster2
tma1so$celltype[suspected_tumor] <- "tumor"
tma1so$celltype[suspected_rbc] <- "erythrocyte"
tma1so$celltype[(tma1so$insitutype_sup_cluster == "Hepatocytes")] <- "tumor"
tma1so$celltype[(tma1so$insitutype_sup_cluster == "Cholangiocytes")] <- "cholangiocyte"

DimPlot(tma1so, group.by = "celltype", raster = F, order = F) + ggprism::scale_color_prism()

SaveSeuratRds(tma1so, file = "hcc_tma1_final.RDS")


# TMA2 -------------------------------------------------------------------------
cidx3 <- !(tma1so$celltype %in% c("tumor", "erythrocyte"))
newref <- InSituType::getRNAprofiles(x = t(tma1so@assays$RNA@counts[,cidx3]), 
                                     neg = t(tma1so@assays$NP@counts[,cidx3]), 
                                     clust = tma1so$celltype[cidx3])

# Identifying tumor
tma2so <- FindNeighbors(object = tma2so, dims = 1:25)
tma2so <- FindClusters(tma2so, algorithm = 1, resolution = seq(0.5, 1.5, 0.1))
DimPlot(tma2so, raster = F, order = T,  label = T, group.by = "RNA_snn_res.0.7") + NoLegend()

FeaturePlot(object = tma2so, features = "Mean.PanCK", order = T, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma2so, features = "Mean.CD45", order = F, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma2so, features = "Mean.CD68_CK8_18", order = T, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma2so, features = "Area.um2", order = T, raster = F)  + 
  scale_color_viridis_c(option = "turbo") +
  theme_void()

FeaturePlot(object = tma2so, features = "APOA1", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma2so, features = "GC", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma2so, features = "VTN", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma2so, features = "ARG1", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma2so, features = "FGG", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()
FeaturePlot(object = tma2so, features = "GPC3", order = T, raster = F) + scale_color_viridis_c(option = "turbo") + theme_void()

tma2so <- AddModuleScore(tma2so, features = list(c("APOA1", "GC", "VTN", "ARG1", "FGG", "GPC3")), nbin = 5, name = "tumor.score")
FeaturePlot(tma2so, features = "tumor.score1", order = T) + scale_color_viridis_c()

suspected_tumor <- (tma2so$RNA_snn_res.0.7 %in% c("1", "2", "3", "13", "10", "11", "6", "8", "4", "7"))
tma2so$suspected_tumor <- NA
tma2so$suspected_tumor[suspected_tumor] <- "tumor"
DimPlot(tma2so, group.by = "suspected_tumor", raster = F, order = T) 

# Identifying Erythrocytes
# Genes with three blue bits in their barcodes should all be expressed in them
tma2so <- AddModuleScore(tma2so, features = list(bluetargets), name = "Autofluorescence.Score", nbin = 5)
FeaturePlot(tma2so, features = "Autofluorescence.Score1", raster = F, order = F) + scale_color_viridis_c()

suspected_rbc <- (tma2so$RNA_snn_res.0.7 == "14")
tma2so$suspected_rbc <- NA
tma2so$suspected_rbc[suspected_rbc] <- "rbc"
DimPlot(tma2so, group.by = "suspected_rbc", raster = F, order = T) 

# Defining cohort data
tma2so$sdimx <- tma2so$CenterX_global_px*(0.51/4256)
tma2so$sdimy <- tma2so$CenterY_global_px*(0.51/4256)

cidx <- (is.na(tma2so$suspected_tumor) & is.na(tma2so$suspected_rbc))
neighbors <- InSituCor:::nearestNeighborGraph(x = tma2so$sdimx, y = tma2so$sdimy, N = 50, subset = tma2so$core) 
neighborexpression <- InSituCor:::get_neighborhood_expression(counts = tma2so@assays$RNA@counts |> t(), neighbors = neighbors) 
neighborhoodPCs <- irlba::prcomp_irlba(neighborexpression, n = 10)$x

cohortdata <- tma2so@meta.data[cidx,] |> 
  dplyr::select(Area, Mean.PanCK, Mean.CD45, Mean.CD68_CK8_18, Mean.DAPI, Mean.CD298_B2M)
cohort <- InSituType::fastCohorting(mat = cbind(cohortdata, neighborhoodPCs[cidx,]))

neg <- tma2so@assays$NP@counts

# Running the function
sup <- insitutype(
  x = tma2so@assays$RNA@counts[,cidx] |> t(),
  neg = neg[,cidx] |> colMeans(),
  assay_type = "RNA", 
  n_clusts=0, 
  anchors = NULL, 
  reference_profiles = meanrefcts,
  reference_sds = sdrefcts,
  cohort = cohort,
  update_reference_profiles = T, 
  rescale = T, 
  refit = F
) 

tma2so$insitutype_sup_cluster <- NA
tma2so$insitutype_sup_cluster[cidx] <- sup$clust

tma2so$insitutype_sup_prob <- NA
tma2so$insitutype_sup_prob[cidx] <- sup$prob

DimPlot(tma2so, group.by = "insitutype_sup_cluster", raster = F, order = F) + ggprism::scale_color_prism()

cidx2 <- !is.na(tma2so$insitutype_sup_cluster)

cohortdata2 <- tma2so@meta.data[cidx2,] |> 
  dplyr::select(Area, Mean.PanCK, Mean.CD45, Mean.CD68_CK8_18, Mean.DAPI, Mean.CD298_B2M)
cohort2 <- InSituType::fastCohorting(mat = cbind(cohortdata2, neighborhoodPCs[cidx2,]))

sup2 <- insitutype(
  x = tma2so@assays$RNA@counts[,cidx2] |> t(),
  neg = neg[,cidx2] |> colMeans(),
  assay_type = "RNA", 
  n_clusts=0, 
  anchors = NULL, 
  reference_profiles = newref,
  cohort = cohort2,
  update_reference_profiles = T, 
  rescale = T, 
  refit = F
) 

tma2so$insitutype_sup_cluster2 <- NA
tma2so$insitutype_sup_cluster2[cidx2] <- sup2$clust

tma2so$insitutype_sup_prob2 <- NA
tma2so$insitutype_sup_prob2[cidx2] <- sup2$prob

DimPlot(tma2so, group.by = "insitutype_sup_cluster2", raster = F, order = F) + ggprism::scale_color_prism()

tma2so$celltype <- tma2so$insitutype_sup_cluster2
tma2so$celltype[suspected_tumor] <- "tumor"
tma2so$celltype[suspected_rbc] <- "erythrocyte"
tma2so$celltype[(tma2so$insitutype_sup_cluster == "Hepatocytes")] <- "tumor"

DimPlot(tma2so, group.by = "celltype", raster = F, order = F) + ggprism::scale_color_prism()

SaveSeuratRds(tma2so, file = "hcc_tma2_final.RDS")


# Merging ----------------------------------------------------------------------
options(Seurat.object.assay.version = "v3")
tmas <- CreateSeuratObject(counts = cbind(tma1so@assays$RNA@counts, tma2so@assays$RNA@counts),
                           meta.data = rbind(tma1so@meta.data, tma2so@meta.data), project = "hcc_tma_project_final")
negs <- cbind(tma1so@assays$NP@counts, tma2so@assays$NP@counts)
scontrs <- cbind(tma1so@assays$FC@counts, tma2so@assays$FC@counts)
tmas[["NP"]] <- CreateAssayObject(counts = negs)
tmas[["FC"]] <- CreateAssayObject(counts = scontrs)
tmas$slide <- paste("hcc_tma", str_split(string = tmas$patient, pattern = "c|p", simplify = T)[,1], sep = "")

tmas <- NormalizeData(tmas, normalization.method = "LogNormalize", scale.factor = 1000, assay = "RNA")
tmas <- FindVariableFeatures(tmas, selection.method = "dispersion", assay = "RNA", nfeatures = nrow(tmas))
tmas <- ScaleData(tmas, assay = "RNA")
tmas <- RunPCA(tmas, npcs = 50, assay = "RNA")
ElbowPlot(tmas, ndims = 50)
tmas <- RunUMAP(tmas, dims = 1:25, assay = "RNA")

DimPlot(tmas, group.by = "celltype", raster = F) + 
  theme_void() + ggprism::scale_color_prism()
DimPlot(tmas, group.by = "slide", raster = F) + 
  theme_void() + 
  NoLegend()

SaveSeuratRds(tmas, "hcc_tmas_final.RDS")

pdf("2-celltyping/celltyping_feature_plots.pdf", width = 8, height = 8)
FeaturePlot(tmas, "CD8A", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "CD3E", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "CD163", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "IGHG1", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "PECAM1", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "COL1A1", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "MMP7", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "ARG1", order = T, raster = F) + scale_color_viridis_c()
FeaturePlot(tmas, "APOA1", order = T, raster = F) + scale_color_viridis_c()
dev.off()

library(lmerTest)
library(presto)

psb <- presto::collapse_counts(counts_mat = tmas@assays$RNA@counts, 
                               meta_data = tmas@meta.data,
                               get_norm = F,
                               varnames = c("celltype", "patient"))
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())

presto_res <- presto.presto(
  y ~ 1 + (1|celltype) + (1|celltype:patient) + (1|patient) + offset(logUMI), 
  psb$meta_data, 
  psb$counts_mat,
  size_varname = 'logUMI', 
  effects_cov = c('celltype'),
  ncore = 10, 
  min_sigma = .05, 
  family = 'poisson',
  nsim = 1000
) 

contrasts_mat <- make_contrast.presto(presto_res, 'celltype')
effects_marginal <- contrasts.presto(presto_res, contrasts_mat, one_tailed = TRUE) %>% 
  dplyr::mutate(cluster = contrast) %>% 
  dplyr::mutate(
    ## convert stats to log2 for interpretability 
    logFC = sign(beta) * log2(exp(abs(beta))),
    SD = log2(exp(sigma)),
    zscore = logFC / SD
  ) %>% 
  dplyr::select(cluster, feature, logFC, SD, zscore, pvalue) %>% 
  arrange(pvalue)
effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, method = 'BH')

effects_marginal |> filter(fdr < 0.05 & logFC > 1) |> arrange(cluster, desc(logFC)) |> 
  ggplot() + geom_point(mapping = aes(x = logFC, y = -log10(fdr)))

pdf("2-celltyping/celltyping_volcano_plots.pdf", height = 40, width = 8)
p <- effects_marginal |> 
  ggplot() + 
  geom_point(mapping = aes(x = logFC, y = -log10(fdr))) + 
  geom_hline(yintercept = -log10(0.05)) + 
  geom_vline(xintercept = 1.5) + 
  ggrepel::geom_text_repel(data = effects_marginal |> filter(fdr < 0.05 & logFC > 1.5), 
                           mapping = aes(x = logFC, y = -log10(fdr), label = feature), 
                           color = "red", size = 2, max.overlaps = 45) + 
  ggthemes::theme_few() + 
  labs(title = "Cell Type Markers") +
  facet_wrap(.~cluster, ncol = 2, scales = "free")
p
dev.off()

write.csv(x = effects_marginal, file = "2-celltyping/celltyping_DE_results.csv")

pdf("2-celltyping/celltyping_dotplot.pdf", height = 10, width = 8)

# Getting average profiles
profiles <- InSituType:::Estep(counts = (tmas@assays$RNA@counts |> t()), 
                               clust = tmas$celltype,
                               neg =  tmas@assays$NP@counts |> colMeans())$profiles

# Getting the proportion of cells that are positive for each gene
meanexpressing <- profiles * NA
for (cell in colnames(meanexpressing)) {
  tempmat <- (tmas@assays$RNA@counts)[rownames(meanexpressing), (tmas$celltype == cell)]
  meanexpressing[, cell] <- Matrix::rowMeans(tempmat > 0)
}

# Genes that act as markers 
markers <- (apply(profiles, 1, max) > 0.1) & 
  (apply(profiles, 1, max) > 2 * apply(profiles, 1, function(x){
    x[order(x, decreasing = T)[2]]}))

# Genes that act as cell type markers
celltypemarkers <- list()
for (cell in colnames(profiles)) {
  celltypemarkers[[cell]] <- names(which(
    (profiles[, cell] > 0.1) & 
      (profiles[, cell] > 1 * apply(profiles, 1, function(x){
        x[order(x, decreasing = T)[2]]}))
  ))
  
  celltypemarkers[[cell]] <- celltypemarkers[[cell]][
    order(profiles[celltypemarkers[[cell]], cell] / 
            apply(profiles[celltypemarkers[[cell]], setdiff(colnames(profiles), cell), drop = F], 1, max),
          decreasing = T)[1:pmin(5, length(celltypemarkers[[cell]]))]]
}

df = data.frame(gene = rep(rownames(profiles), ncol(profiles)),
                clust = rep(colnames(profiles), each = nrow(profiles)),
                mean = as.vector(profiles),
                prop.expressing = as.vector(meanexpressing))
df <- df[is.element(df$gene, unlist(celltypemarkers)), ]
df$scaled.mean = df$mean
for (gene in unique(df$gene)) {
  df$scaled.mean[df$gene == gene] = df$scaled.mean[df$gene == gene] / max(df$scaled.mean[df$gene == gene])
}

geneorder <- unlist(celltypemarkers) |> na.omit()
df$gene = factor(df$gene, levels = geneorder)

cellorder <- names(celltypemarkers)
df$clust = factor(df$clust, levels = cellorder)

ggplot(df, aes(x=clust, y = gene, color = scaled.mean, size = prop.expressing)) + 
  geom_point() + scale_color_viridis_c(option = "A", direction = -1) + ggthemes::theme_few() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y = element_text(size = 7))

dev.off()

