rm(list = ls())
.libPaths()

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(magrittr)
library(ComplexHeatmap)

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
                   "Treg"="pink3",
                   "mDC"="yellow3",
                   "T CD8 naive"="orange",
                   "T CD4 memory"="brown2",
                   "monocyte"="darkblue",     
                   "T CD8 memory"="beige",
                   "NK"="tan",
                   "pDC"="orange3", 
                   "T CD4 naive"="turquoise", 
                   "erythrocyte"="pink")  
tmas$core <- paste(substr(x = tmas$slide, start = 8, stop = 8), tmas$core, sep = "")
tmas$fov <- paste(substr(x = tmas$slide, start = 8, stop = 8), tmas$fov, sep = "f")
tmas$celltype %<>% as.factor()
tmas$patient %<>% as.character() %<>% as.factor()

## FIGURE 1 --------------------------------------------------------------------
# D
# tma2fov60 <- jpeg::readJPEG("7-figures/fig1d/fig1d_left.jpg")
txoi <- arrow::open_csv_dataset("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc_tma2/flatFiles/hcc_tma2/hcc_tma2_tx-hive") |> 
  filter(fov == 60) |> filter(target %in% c("GC", "APOC1", "TTR", "CD4", "HLA-DRA", "HLA-DRB", "CD68", "C1QC", "SPOCK2", "PTGDS", "IGKC")) |> 
  collect()
txoi$tx_category <- case_when( 
                              (txoi$target %in% c("GC", "APOC1", "TTR")) ~ "Tumor", 
                              (txoi$target) %in% c("CD4", "HLA-DRA", "HLA-DRB", "CD68", "C1QC", "SPOCK2", "PTGDS", "IGKC") ~ "Immune"
                              )
# fig1d_topright <- ggplot(mapping = aes(x = x_global_px, y = y_global_px, color = tx_category)) + 
#   ggpubr::background_image(tma2fov60) + 
#   geom_point(data = txoi |> filter(tx_category == "Tumor"), size = 0.1) + 
#   annotate(geom = "segment", x = 48500, xend = 48500+100*(4256/1000/0.51), y = 3500, yend = 3500, color = "white") +
#   annotate(geom = "text", x = 48800, y = 3650, color = "white", label = "100 µm", size = 4, fontface = "bold") +
#   scale_color_manual(values = c("Immune"="dodgerblue", "Tumor"="yellow3")) +
#   coord_fixed(expand = F) + 
#   theme_void() + 
#   theme(strip.text = element_blank()) + 
#   NoLegend()
# 
# fig1d_bottomright <- ggplot(mapping = aes(x = x_global_px, y = y_global_px, color = tx_category)) + 
#   ggpubr::background_image(tma2fov60) + 
#   geom_point(data = txoi |> filter(tx_category == "Immune"), size = 0.1) + 
#   annotate(geom = "segment", x = 48500, xend = 48500+100*(4256/1000/0.51), y = 3500, yend = 3500, color = "white") +
#   annotate(geom = "text", x = 48800, y = 3650, color = "white", label = "100 µm", size = 4, fontface = "bold") +
#   scale_color_manual(values = c("Immune"="dodgerblue", "Tumor"="yellow3")) +
#   coord_fixed(expand = F) + 
#   theme_void() + 
#   theme(strip.text = element_blank()) + 
#   NoLegend()
# 
# fig1d_topright
# #ggsave(filename = "7-figures/fig1d/fig1d_topright.pdf", device = "pdf", width = 6, height = 6, dpi = "screen")
# fig1d_bottomright
# #ggsave(filename = "7-figures/fig1d/fig1d_bottomright.pdf", device = "pdf", width = 6, height = 6, dpi = "screen")

tma2fov60 <- jpeg::readJPEG("7-figures/fig1d/fig1d_left_4X_downsampled.jpg")
fig1d_topright <- ggplot(mapping = aes(x = x_global_px, y = y_global_px, color = tx_category)) + 
  ggpubr::background_image(tma2fov60) + 
  geom_point(data = txoi |> filter(tx_category == "Tumor"), size = 0.1) + 
  annotate(geom = "segment", x = 48500, xend = 48500+100*(4256/1000/0.51), y = 3500, yend = 3500, color = "white") +
  annotate(geom = "text", x = 48800, y = 3650, color = "white", label = "100 µm", size = 4, fontface = "bold") +
  scale_color_manual(values = c("Immune"="dodgerblue", "Tumor"="yellow3")) +
  coord_fixed(expand = F) + 
  theme_void() + 
  theme(strip.text = element_blank()) + 
  NoLegend()

fig1d_bottomright <- ggplot(mapping = aes(x = x_global_px, y = y_global_px, color = tx_category)) + 
  ggpubr::background_image(tma2fov60) + 
  geom_point(data = txoi |> filter(tx_category == "Immune"), size = 0.1) + 
  annotate(geom = "segment", x = 48500, xend = 48500+100*(4256/1000/0.51), y = 3500, yend = 3500, color = "white") +
  annotate(geom = "text", x = 48800, y = 3650, color = "white", label = "100 µm", size = 4, fontface = "bold") +
  scale_color_manual(values = c("Immune"="dodgerblue", "Tumor"="yellow3")) +
  coord_fixed(expand = F) + 
  theme_void() + 
  theme(strip.text = element_blank()) + 
  NoLegend()

fig1d_topright
#ggsave(filename = "7-figures/fig1d/fig1d_topright_4Xdownsampled.pdf", device = "pdf", width = 6, height = 6, dpi = "screen")
fig1d_bottomright
#ggsave(filename = "7-figures/fig1d/fig1d_bottomright_4Xdownsampled.pdf", device = "pdf", width = 6, height = 6, dpi = "screen")

## FIGURE 2 --------------------------------------------------------------------
# A
tma1 <- LoadSeuratRds("hcc_tma1_final.RDS")
tma2 <- LoadSeuratRds("hcc_tma2_final.RDS")
fig2a_1 <- 
(FeaturePlot(object = tma1, features = c("tumor.score1"), raster = T) + 
  scale_color_viridis_c() + 
   scale_x_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
   scale_y_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
   coord_fixed() +
   theme_classic() +
   theme(
     axis.text = element_blank(), 
     axis.title.x = element_text(face = "bold", size = 12),
     axis.title.y = element_text(face = "bold", size = 12),
     axis.ticks = element_blank(),
     axis.title = element_text(hjust = 0.125),
     axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
     legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(10, units = "pt"), legend.position = "inside", legend.position.inside = c(0.8, 0.8),
     plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
     ) +
  labs(title = "Tumor\nMetascore")) # /
  # (FeaturePlot(object = tma2, features = c("tumor.score1"), raster = T) + 
  #    scale_color_viridis_c() + 
  #    scale_x_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  #    scale_y_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  #    coord_fixed() +
  #    theme_classic() +
  #    theme(
  #      axis.text = element_blank(), 
  #      axis.title.x = element_text(face = "bold", size = 10),
  #      axis.title.y = element_text(face = "bold", size = 10),
  #      axis.ticks = element_blank(),
  #      axis.title = element_text(hjust = 0.125),
  #      axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
  #      legend.text = element_text(face = "bold", size = 5), legend.key.size = unit(5, units = "pt"),
  #      plot.title = element_blank()
  #    )
  #  )
fig2a_1

fig2a_2 <- 
  (FeaturePlot(object = tma1, features = c("Autofluorescence.Score1"), raster = T) + 
     scale_color_viridis_c() + 
     scale_x_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
     scale_y_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
     coord_fixed() +
     theme_classic() +
     theme(
       axis.text = element_blank(), 
       axis.title.x = element_text(face = "bold", size = 12),
       axis.title.y = element_text(face = "bold", size = 12),
       axis.ticks = element_blank(),
       axis.title = element_text(hjust = 0.125),
       axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
       legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(10, units = "pt"), legend.position = "inside", legend.position.inside = c(0.8, 0.8),
       plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
     ) +
     labs(title = "Autofluorescence\nMetascore")) # /
  # (FeaturePlot(object = tma2, features = c("Autofluorescence.Score1"), raster = T) + 
  #    scale_color_viridis_c() + 
  #    scale_x_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  #    scale_y_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  #    coord_fixed() +
  #    theme_classic() +
  #    theme(
  #      axis.text = element_blank(), 
  #      axis.title.x = element_text(face = "bold", size = 10),
  #      axis.title.y = element_text(face = "bold", size = 10),
  #      axis.ticks = element_blank(),
  #      axis.title = element_text(hjust = 0.125),
  #      axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
  #      legend.text = element_text(face = "bold", size = 5), legend.key.size = unit(5, units = "pt"),
  #      plot.title = element_blank()
  #    )
  #  )
fig2a_2

fig2a_3 <- 
  (DimPlot(object = tma1, group.by = "celltype", raster = T) + 
     scale_color_manual(values = c("tumor"="green4", "erythrocyte"="pink")) + 
     scale_x_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
     scale_y_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
     coord_fixed() +
     theme_classic() +
     theme(
       axis.text = element_blank(), 
       axis.title.x = element_text(face = "bold", size = 12),
       axis.title.y = element_text(face = "bold", size = 12),
       axis.ticks = element_blank(),
       axis.title = element_text(hjust = 0.125),
       axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
       legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(10, units = "pt"), legend.position = "inside", legend.position.inside = c(0.8, 0.7),
       plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
     ) +
     labs(title = "Tumor Cells &\nErythrocytes")) # /
  # (DimPlot(object = tma2, group.by = "celltype", raster = T) + 
  #    scale_color_manual(values = c("tumor"="green4", "erythrocyte"="pink")) + 
  #    scale_x_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  #    scale_y_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  #    coord_fixed() +
  #    theme_classic() +
  #    theme(
  #      axis.text = element_blank(), 
  #      axis.title.x = element_text(face = "bold", size = 10),
  #      axis.title.y = element_text(face = "bold", size = 10),
  #      axis.ticks = element_blank(),
  #      axis.title = element_text(hjust = 0.125),
  #      axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
  #      legend.text = element_text(face = "bold", size = 5), legend.key.size = unit(5, units = "pt"), legend.position = "top",
  #      plot.title = element_blank()
  #    ) + 
  #    NoLegend()
  # )
fig2a_3

fig2a_4 <- 
  (DimPlot(object = tma1, group.by = "celltype", raster = T) + 
     scale_color_manual(values = celltype_cols) + 
     scale_x_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
     scale_y_continuous(breaks = max(tma1@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
     coord_fixed() +
     theme_classic() +
     theme(
       axis.text = element_blank(), 
       axis.title.x = element_text(face = "bold", size = 12),
       axis.title.y = element_text(face = "bold", size = 12),
       axis.ticks = element_blank(),
       axis.title = element_text(hjust = 0.125),
       axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
       legend.text = element_text(face = "bold", size = 5), legend.key.size = unit(5, units = "pt"),
       plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
     ) +
     labs(title = "Final\nCell Types")) # /
  # (DimPlot(object = tma2, group.by = "celltype", raster = T) + 
  #    scale_color_manual(values = celltype_cols) + 
  #    scale_x_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  #    scale_y_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  #    coord_fixed() +
  #    theme_classic() +
  #    theme(
  #      axis.text = element_blank(), 
  #      axis.title.x = element_text(face = "bold", size = 10),
  #      axis.title.y = element_text(face = "bold", size = 10),
  #      axis.ticks = element_blank(),
  #      axis.title = element_text(hjust = 0.125),
  #      axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
  #      legend.text = element_text(face = "bold", size = 5), legend.key.size = unit(5, units = "pt"),
  #      plot.title = element_blank()
  #    )
  # )
fig2a_4

(fig2a_1|(patchwork::plot_spacer())|fig2a_2)
#ggsave(filename = "7-figures/fig2a_topleft.pdf", device = "pdf", width = 8, height = 4.5)

fig2a_3
#ggsave(filename = "7-figures/fig2a_topright.pdf", device = "pdf", width = 4, height = 4.5)

fig2a_4 + NoLegend()
#ggsave(filename = "7-figures/fig2a_bottomright.pdf", device = "pdf", width = 4, height = 4.5)

# B
fig2b_left <- DimPlot(object = tmas, group.by = "celltype", raster = T) + 
  scale_color_manual(values = celltype_cols) + 
  scale_x_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tma2@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(5, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 15), legend.key.size = unit(15, units = "pt"),
    plot.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
fig2b_left
#ggsave(filename = "7-figures/fig2b_left.pdf", device = "pdf", width = 8, height = 8)

fig2b_right <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "2f60"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype)) + 
  scale_color_manual(values = celltype_cols) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
        ) + 
  labs(title = "TMA2 FOV60") +
  NoLegend()
fig2b_right
#ggsave(filename = "7-figures/fig2b_right.pdf", device = "pdf", width = 6, height = 6)

# C 
tmas$celltype %<>% as.character()
profiles <- InSituType:::Estep(counts = (tmas@assays$RNA@counts |> t()), 
                               clust = tmas$celltype,
                               neg =  tmas@assays$NP@counts |> colMeans())$profiles

meanexpressing <- profiles * NA
for (cell in colnames(meanexpressing)) {
  tempmat <- (tmas@assays$RNA@counts)[rownames(meanexpressing), (tmas$celltype == cell)]
  meanexpressing[, cell] <- Matrix::rowMeans(tempmat > 0)
}

markers <- (apply(profiles, 1, max) > 0.1) & 
  (apply(profiles, 1, max) > 2 * apply(profiles, 1, function(x){
    x[order(x, decreasing = T)[2]]}))

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

fig2c <- ggplot(df, aes(x=gene, y = clust, fill = scaled.mean, size = prop.expressing)) + 
  geom_point(shape = 21, color = "black") + scale_fill_viridis_c(option = "A", direction = -1) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 8), axis.text.y = element_text(hjust = 1), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  labs(size = "Proportion\nPositive", fill = "Scaled Mean\nExpression")
fig2c
#ggsave(filename = "7-figures/fig2c.pdf", device = "pdf", width = 14, height = 6)

# D
cellgroup_map <- c(
  "B-cell"="lymphoid",
  "cholangiocyte"="cholangiocyte", 
  "endothelial"="stromal", 
  "erythrocyte"="erythrocyte", 
  "fibroblast"="stromal",
  "macrophage"="myeloid",
  "mast"="myeloid",
  "mDC"="myeloid",
  "monocyte"="myeloid",
  "neutrophil"="myeloid",
  "NK"="lymphoid",
  "pDC"="lymphoid",
  "plasmablast"="lymphoid",
  "T CD4 memory"="lymphoid",
  "T CD4 naive"="lymphoid",
  "T CD8 memory"="lymphoid",
  "T CD8 naive"="lymphoid",
  "Treg"="lymphoid",
  "tumor"="tumor"
)

tmas$cellgroup <- plyr::mapvalues(x = tmas$celltype, from = names(cellgroup_map), to = cellgroup_map)
tmas$cellgroup <- factor(tmas$cellgroup, levels = c("tumor", "stromal", "erythrocyte", "cholangiocyte", "lymphoid", "myeloid"))
cellcts <- tmas@meta.data |>
  dplyr::group_by(patient, fov, cellgroup, .drop = F) |>
  tally() |> dplyr::group_by(patient, fov, .drop = F) |> dplyr::mutate(total = sum(n), prop = n/sum(n)) |> 
  as.data.frame()

fig2d <- ggplot(cellcts) + 
  geom_bar(mapping = aes(x = fov, y = prop, fill = cellgroup, group = fov), stat = "identity", color = NA) + 
  scale_fill_manual(values = c("tumor" = "green4", "stromal" = "grey", "erythrocyte" = "pink", "cholangiocyte" = "dodgerblue", "lymphoid" = "violet", "myeloid" = "orange")) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7), 
        legend.title = element_blank(), legend.text = element_text(face = "bold", size = 20), legend.position = "top") + 
  labs(x = "FOV", y = "proportion") + 
  guides(fill = guide_legend(nrow = 1))
fig2d
#ggsave(filename = "7-figures/fig2d.pdf", device = "pdf", width = 15, height = 5)

## FIGURE 3 --------------------------------------------------------------------
# A
ptkey <- read.csv("3-DE-analysis/cosmx-patient-data-de-identified.csv", row.names = 1)
ptkey$patient_deid <- factor(ptkey$patient_deid, levels = dplyr::arrange(ptkey, dplyr::desc(line1orf1_tumor_cpm)) |> pull(patient_deid))
fig3a <- ggplot() + 
  geom_bar(data = ptkey, mapping = aes(x = patient_deid, y = line1orf1_tumor_cpm, fill = line1orf1_group), stat = "identity") + 
  geom_hline(yintercept = quantile(x = ptkey$line1orf1_tumor_cpm, probs = 2/3), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("high" = "firebrick", "low" = "dodgerblue")) + 
  ggthemes::theme_par() +
  ggpubr::labs_pubr() + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.text = element_text(face = "bold"), legend.position = "inside", legend.position.inside = c(0.85, 0.85), legend.title = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8), 
        axis.text.y = element_text(hjust = 1)) + 
  labs(y = "LINE1-ORF1 CPM in Tumor Cells", title = "LINE1-ORF1 High and Low Patients")
fig3a
#ggsave(filename = "7-figures/fig3a.pdf", device = "pdf", width = 6, height = 6)

# fig3a <- ggplot() + 
#   stat_summary(data = ptkey, mapping = aes(x = line1orf1_group, y = line1orf1_tumor_cpm), geom = "crossbar", fun = "mean") +
#   geom_jitter(data = ptkey, mapping = aes(x = line1orf1_group, y = line1orf1_tumor_cpm, fill = line1orf1_group, shape = line1orf1_group), height = 0, width = 0.25, size = 3) + 
#   geom_hline(yintercept = quantile(x = ptkey$line1orf1_tumor_cpm, probs = 2/3), linetype = "dashed", color = "black") + 
#   scale_fill_manual(values = c("high" = "firebrick", "low" = "dodgerblue")) + 
#   scale_shape_manual(values = c("high"=21, "low"=22)) +
#   ggthemes::theme_par() +
#   ggpubr::labs_pubr() + 
#   scale_y_continuous(limits = c(0, 35000), breaks = c(0, 10000, 20000, 30000)) + 
#   theme(legend.text = element_text(face = "bold"), legend.position = "inside", legend.position.inside = c(0.8, 0.9), legend.title = element_blank(), 
#         axis.title.x = element_blank(), 
#         axis.text.y = element_text(hjust = 1)) + 
#   labs(y = "LINE1-ORF1 CPM in Tumor Cells", title = "LINE1-ORF1 High and Low Patients")
# fig3a

# B
tumorres <- read.csv("3-DE-analysis/line1_patient_level_de_analysis_tumor_results.csv", row.names = 1)
fig3b <- ggplot() + 
  # Lines ######################################################################
  geom_hline(yintercept = -log10(0.05), linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  # NS points ##################################################################
  geom_point(data = tumorres |> filter(p.adj > 0.05), 
             mapping = aes(x = logFC, y = -log10(p.adj)), color = "grey") +
  # Left points ################################################################ 
  geom_point(data = tumorres |> filter(p.adj < 0.05 & logFC < 0), 
             mapping = aes(x = logFC, y = -log10(p.adj)), color = "dodgerblue") +
  ggrepel::geom_text_repel(data = filter(tumorres, (p.adj < 0.05) & (logFC < 0)), 
                           mapping = aes(x = logFC, y = -log10(p.adj), label = target, fontface = "bold"), 
                           size = 3, max.overlaps = 15, segment.size = 0.25, force_pull = 1, min.segment.length = 0, ylim = c(1.5, 10), xlim = c(0, -1.5), box.padding = 0.5) +
  # Right points ###############################################################
  geom_point(data = tumorres |> filter(p.adj < 0.05 & logFC > 0), 
             mapping = aes(x = logFC, y = -log10(p.adj)), color = "firebrick") +
  ggrepel::geom_text_repel(data = filter(tumorres, (p.adj < 0.05) & (logFC > 0) & target %in% c("LINE1-ORF1", "LINE1-ORF2", "HERVK", "HSATII", "POU5F1", "LEFTY1", "HCAR2/3", "EFNA4", "EFNA5", "TNF", "CYTOR", "HEY1", "CLDN4", "TGFB3", "TNFRSF18", "TIGIT", "PDCD1", "TWIST2")), 
                           mapping = aes(x = logFC, y = -log10(p.adj), label = target, fontface = "bold"), 
                           size = 3, max.overlaps = 15, segment.size = 0.25, force_pull = 1, min.segment.length = 0, ylim = c(1.5, 10), xlim = c(0, 2), box.padding = 0.5) +
  # Arrows #####################################################################
  annotate(geom = "segment", x = 0, xend = 0.5, y = -0.5, yend = -0.5, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -0.5, label = "Up in High", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -0.5, yend = -0.5, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = -0.6, y = -0.5, label = "Up in Low", hjust = 1, fontface = "bold") +
  # Theme ######################################################################
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  scale_y_continuous(limits = c(-1, 10), breaks = 0:10) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Tumor Cell DEGs: LINE1-ORF1 High vs Low Patients")

fig3b
#ggsave(filename = "7-figures/fig3b.pdf", device = "pdf", width = 8, height = 8)

# C
tmp_norm <- tmas@assays$RNA@data[,tmas$celltype == "tumor"]
tmp_meta <- tmas@meta.data |> filter(celltype == "tumor")
tmp_means_mm <- model.matrix(~0+patient, data = tmp_meta)
tmp_means <- (tmp_norm %*% tmp_means_mm)
tmp_means %<>% as.matrix()
tmp_means <- sweep(x = tmp_means, MARGIN = 2, STATS = colSums(tmp_means_mm), FUN = "/")
colnames(tmp_means) <- gsub(pattern = "patient", replacement = "", x = colnames(tmp_means))
tmp_means <- apply(X = tmp_means, MARGIN = 1, FUN = scale) |> t()
colnames(tmp_means) <- colnames(tmp_means) <- gsub(pattern = "patient", replacement = "", x = colnames(tmp_means_mm))

#pdf(file = "7-figures/fig3c.pdf", width = 8, height = 8)
mat <- tmp_means[c("LINE1-ORF1", "LINE1-ORF2", "HERVK", "HSATII", "POU5F1", "LEFTY1", "HCAR2/3", "EFNA4", "EFNA5", "TNF", "CYTOR", "HEY1", "CLDN4", "TGFB3", "TNFRSF18", "TIGIT", "PDCD1", "TWIST2", tumorres[tumorres$p.adj < 0.05 & tumorres$logFC < 0,]$target), levels(ptkey$patient_deid)]
Heatmap(matrix = mat, 
        cluster_columns = F, 
        column_split = (ptkey |> arrange(line1orf1_group) |> pull(line1orf1_group)), 
        column_title_gp = gpar(fontface = "bold"),
        name = "Mean Expression", 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(3, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
        column_names_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(group = ptkey |> arrange(line1orf1_group) |> pull(line1orf1_group), 
                                           col = list(group = c("low"="dodgerblue", "high"="firebrick")), 
                                           show_legend = F, 
                                           annotation_name_gp = gpar(col = NA)),
        #col = circlize::colorRamp2(breaks = c(-4, -2, 0, 2, 4), colors = c("blue", "blue", "white", "red", "red")), 
        heatmap_legend_param = list(direction = "horizontal")
) |> 
  draw(heatmap_legend_side = "bottom")
#dev.off()

# D
tumor_cells <- (tmas$celltype == "tumor")
tmas$line1orf1 <- tmas@assays$RNA@data["LINE1-ORF1",]
meta <- tmas@meta.data[tumor_cells,]
meta$patient <- as.character(meta$patient) |> as.factor()
meta <- dplyr::group_by(meta, patient) |> dplyr::mutate(cut = quantile(line1orf1, 2/3))
meta$line1orf1_status <- ifelse(test = (meta$line1orf1 < meta$cut), yes = "L1low", no = "L1high")
meta$patient <- factor(x = meta$patient, levels = dplyr::arrange(ptkey, dplyr::desc(line1orf1_tumor_cpm)) |> pull(patient_deid))
fig3d <- ggplot() + 
  scattermore::geom_scattermore(data = meta, mapping = aes(x = patient, y = line1orf1), pointsize = 0, position = "jitter", alpha  = 0.25) +
  geom_boxplot(data = meta, mapping = aes(x = patient, y = line1orf1, fill = line1orf1_status), color = "black", outliers = F, position = "dodge2") +
  stat_summary(data = meta, mapping = aes(x = patient, y = line1orf1), geom = "crossbar", fun = quantile, fun.args = list(2/3)) +
  scale_fill_manual(values = c("L1high" = "firebrick", "L1low" = "dodgerblue")) + 
  ggthemes::theme_par() +
  ggpubr::labs_pubr() + 
  theme(legend.text = element_text(face = "bold", size = 12), legend.position = "inside", legend.position.inside = c(0.9, 0.85), legend.title = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12), 
        axis.text.y = element_text(hjust = 1), plot.title = element_text(size = 20)) + 
  labs(y = "LINE1-ORF1 log2(CPT+1)", title = "LINE1-ORF1 High and Low Tumor Cells")
fig3d
#ggsave(filename = "7-figures/fig3d.pdf", device = "pdf", width = 12, height = 6)

# E
tumorcellres <- read.csv("3-DE-analysis/line1_cell_level_de_analysis_tumor_results.csv", row.names = 1)
fig3e <- ggplot() + 
  # Lines ######################################################################
  geom_hline(yintercept = -log10(0.05), linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  # NS points ##################################################################
  geom_point(data = tumorcellres |> filter(p.adj > 0.05), 
            mapping = aes(x = logFC, y = -log10(p.adj)), color = "grey") +
  # Left points ################################################################ 
  geom_point(data = tumorcellres |> filter(p.adj < 0.05 & logFC < 0), 
             mapping = aes(x = logFC, y = -log10(p.adj)), color = "dodgerblue") +
  ggrepel::geom_text_repel(data = filter(tumorcellres, (p.adj < 0.05) & (logFC < 0)), 
                           mapping = aes(x = logFC, y = -log10(p.adj), label = target, fontface = "bold"), 
                           size = 3, max.overlaps = 55, segment.size = 0.25, force_pull = 1, min.segment.length = 0, box.padding = 0.85, ylim = c(2, 15), xlim = c(-0.05, -1)) +
  # Right points ###############################################################
  geom_point(data = tumorcellres |> filter(p.adj < 0.05 & logFC > 0), 
             mapping = aes(x = logFC, y = -log10(p.adj)), color = "firebrick") +
  ggrepel::geom_text_repel(data = filter(tumorcellres, (p.adj < 0.05) & (logFC > 0)), 
                           mapping = aes(x = logFC, y = -log10(p.adj), label = target, fontface = "bold"), 
                           size = 3, max.overlaps = 55, segment.size = 0.25, force_pull = 1, min.segment.length = 0, box.padding = 1, ylim = c(2, 15), xlim = c(0.05, 1)) +
  # Arrows #####################################################################
  annotate(geom = "segment", x = 0, xend = 0.5, y = -0.5, yend = -0.5, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -0.5, label = "Up in L1high", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -0.5, yend = -0.5, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = -0.6, y = -0.5, label = "Up in L1low", hjust = 1, fontface = "bold") +
  # Theme ######################################################################
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  scale_y_continuous(limits = c(-1, 15), breaks = 0:15) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Tumor Cell DEGs: LINE1-ORF1 High vs Low Cells")
fig3e
#ggsave(filename = "7-figures/fig3e.pdf", device = "pdf", width = 8, height = 8)

# F
library(data.table)
library(DESeq2)
psb <- presto::collapse_counts(counts_mat = tmas@assays$RNA@counts[,tumor_cells], 
                               meta_data = meta, 
                               varnames = c("patient", "line1orf1_status"), 
                               get_norm  = F)
dds <- DESeqDataSetFromMatrix(countData = psb$counts_mat, colData = psb$meta_data, design = ~line1orf1_status+patient)
nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "upperquartile")
sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
sfs <- sfs / exp(mean(log(sfs)))
dds@colData$sizeFactor <- sfs
dd <- DESeq2::counts(object = dds, normalized = T)[c("HSATII", "HERVK", "POU5F1", "LEFTY1", "SERPINA3"),] |> t() |> cbind(psb$meta_data)
dd <- tidyr::pivot_longer(data = dd, cols = c("HSATII", "HERVK", "POU5F1", "LEFTY1", "SERPINA3"), names_to = "gene", values_to = "norm")
dd$gene <- factor(dd$gene, levels = c("HSATII", "HERVK", "POU5F1", "LEFTY1", "SERPINA3"))
fig3f <- ggplot() + 
  geom_line(data = dd, mapping = aes(x = line1orf1_status, y = norm, group = patient), color = "grey") +
  geom_point(data = dd, mapping = aes(x = line1orf1_status, y = norm, group = patient)) + 
  stat_summary(data = dd, mapping = aes(x = line1orf1_status, y = norm, color = line1orf1_status), geom = "crossbar", fun = mean) +
  scale_color_manual(values = c("L1high" = "firebrick", "L1low" = "dodgerblue")) + 
  ggthemes::theme_par() +
  ggpubr::labs_pubr() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(hjust = 1), axis.title.y = element_text(size = 20), 
        strip.text = element_text(face = "bold", size = 20)) + 
  labs(y = "Normalized Expression") + 
  facet_wrap(.~gene, scales = "free_y", nrow = 1) + 
  NoLegend()
fig3f
#ggsave(filename = "7-figures/fig3f.pdf", device = "pdf", width = 11.5, height = 5.5)

# G
library(msigdbr)
c8genesets <- msigdbr(species = "Homo sapiens", category = "C8") |> 
  dplyr::select(gs_name, gene_symbol)
c8genesets %<>% filter(grepl("LIVER", gs_name))
mfgenesets <- msigdbr(species = "Homo sapiens", category = "C5") |> 
  filter(gs_subcat %in% c("GO:MF")) |>
  dplyr::select(gs_name, gene_symbol)

tumorcellres$signed_stat <- tumorcellres$stat*sign(tumorcellres$logFC)
tumgenelist <- dplyr::arrange(tumorcellres |> filter(!(target %in% c("LINE1-ORF1", "LINE1-ORF2", "HSATII", "HERVK", "HERVH"))), dplyr::desc(signed_stat)) |> pull(signed_stat)
names(tumgenelist) <- dplyr::arrange(tumorcellres |> filter(!(target %in% c("LINE1-ORF1", "LINE1-ORF2", "HSATII", "HERVK", "HERVH"))), dplyr::desc(signed_stat)) |> pull(target)

set.seed(2001)
tum_c8 <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = c8genesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_c8@result |> write.csv("3-DE-analysis/tumor_L1high_vs_L1low_c8_significant_gsea_results.csv")
fig3g_topleft <- enrichplot::gseaplot(x = tum_c8, geneSetID = "AIZARANI_LIVER_C14_HEPATOCYTES_2", 
                     color.line = "dodgerblue", color.vline = "black", 
                     title = "AIZARANI_LIVER_C14_HEPATOCYTES_2", by = "runningScore") + 
  theme(plot.title = element_text(size = 20, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 15, face = "bold"))
fig3g_topright <- enrichplot::gseaplot(x = tum_c8, geneSetID = "DESCARTES_FETAL_LIVER_HEPATOBLASTS", 
                                      color.line = "dodgerblue", color.vline = "black", 
                                      title = "DESCARTES_FETAL_LIVER_HEPATOBLASTS", by = "runningScore")

tum_mf <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = mfgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_mf@result |> write.csv("3-DE-analysis/tumor_L1high_vs_L1low_gomf_significant_gsea_results.csv")
fig3g_bottomleft <- enrichplot::gseaplot(x = tum_mf, geneSetID = "GOMF_CELL_ADHESION_MOLECULE_BINDING", 
                     color.line = "dodgerblue", color.vline = "black", 
                     title = "GOMF_CELL_ADHESION_MOLECULE_BINDING", by = "runningScore")
fig3g_bottomright <- enrichplot::gseaplot(x = tum_mf, geneSetID = "GOMF_GROWTH_FACTOR_ACTIVITY", 
                     color.line = "firebrick", color.vline = "black", 
                     title = "GOMF_GROWTH_FACTOR_ACTIVITY", by = "runningScore") + 
  theme(plot.title = element_text(size = 20, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"), axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 15, face = "bold"))

# fig3g <- (fig3g_topleft|fig3g_topright)/(fig3g_bottomleft|fig3g_bottomright)
fig3g <- (fig3g_bottomright)/(fig3g_topleft) # DT wanted to simplify the plot
fig3g
#ggsave("7-figures/fig3g.pdf", device = "pdf", width = 8, height = 8)
  

## FIGURE 4 --------------------------------------------------------------------
# C
neres <- readRDS("6-spatial-analysis/line1_patient_level_ne_analysis_results_final_testing.RDS")
highmat <- neres$obs$high
rownames(highmat) <- gsub(pattern = "celltype", replacement = "", x = rownames(highmat))
lowmat <- neres$obs$low
rownames(lowmat) <- gsub(pattern = "celltype", replacement = "", x = rownames(lowmat))
diffmat <- (neres$obs$high-neres$obs$low)
rownames(diffmat) <- gsub(pattern = "celltype", replacement = "", x = rownames(diffmat))
pvalmat <- neres$pval
rownames(pvalmat) <- gsub(pattern = "celltype", replacement = "", x = rownames(pvalmat))

o <- c("T CD4 naive", "mDC", "Treg", "mast", "NK", "T CD8 naive", "pDC", "T CD4 memory", "B-cell", "monocyte", "T CD8 memory", "neutrophil", "plasmablast", "macrophage", "erythrocyte", "endothelial", "fibroblast", "cholangiocyte", "tumor")
ra <- rowAnnotation(celltype = o, col = list("celltype"=celltype_cols), gp = gpar(col = "white", lwd = 0.5), show_legend = F, annotation_name_gp = gpar(col = NA))
ta <- HeatmapAnnotation(celltype = o, col = list("celltype"=celltype_cols), gp = gpar(col = "white", lwd = 0.5), show_legend = F, annotation_name_gp = gpar(col = NA))
hh <- Heatmap(matrix = highmat[o,o], cluster_rows = F, cluster_columns = F, top_annotation = ta, left_annotation = ra,
              col = circlize::colorRamp2(breaks = c(-5, -2.5, 0, 2.5, 5), colors = c(scales::muted("blue"), "blue", "white", "red", scales::muted("red"))),
        width = ncol(highmat)*unit(5, "mm"),
        height = nrow(highmat)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        column_title = "LINE1-ORF1 High Patients", column_title_gp = gpar(fontface = "bold", col = "firebrick", just = 0.5),
        name = "log2ER_high",
        row_names_side = "left", show_column_names = F, border = T, border_gp = gpar(col = "black"), row_names_gp = gpar(fontface = "bold"),
        heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop", title = "log2ER")
)

lh <- Heatmap(matrix = lowmat[o,o], cluster_rows = F, cluster_columns = F, top_annotation = ta,
        col = circlize::colorRamp2(breaks = c(-5, -2.5, 0, 2.5, 5), colors = c(scales::muted("blue"), "blue", "white", "red", scales::muted("red"))),
        width = ncol(lowmat)*unit(5, "mm"),
        height = nrow(lowmat)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        column_title = "LINE1-ORF1 Low Patients", column_title_gp = gpar(fontface = "bold", col = "dodgerblue", just = 0.5),
        name = "log2ER_low",
        show_row_names = F, show_column_names = F, border = T, border_gp = gpar(col = "black"),
        show_heatmap_legend = F
)

dh <- Heatmap(matrix = diffmat[o,o], cluster_rows = F, cluster_columns = F, top_annotation = ta, left_annotation = ra,
              col = circlize::colorRamp2(breaks = c(-10, -5, 0, 5, 10), colors = c(scales::muted("dodgerblue"), "dodgerblue", "white", "firebrick", scales::muted("firebrick"))),
              width = ncol(diffmat)*unit(5, "mm"),
              height = nrow(diffmat)*unit(5, "mm"), 
              rect_gp = gpar(col = "white", lwd = 0.5), 
              column_title = "High versus Low", column_title_gp = gpar(fontface = "bold", col = "black", just = 0.5),
              name = "log2FC",
              row_names_side = "left", show_column_names = F, border = T, border_gp = gpar(col = "black"), row_names_gp = gpar(fontface = "bold"),
              heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"), 
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(pvalmat[o,o][i, j] < 0.05) {
                  if(diffmat[o,o][i,j] < -3) {
                    grid.text("*", x, y, gp = gpar(col = "white"))
                  }
                  else {
                    grid.text("*", x, y, gp = gpar(col = "black"))
                  }
                } 
                else {
                  grid.text("", x, y)
                }
              }
)

#pdf(file = "7-figures/fig4c.pdf", width = 14, height = 5)
draw(hh + lh, heatmap_legend_side = "bottom")
decorate_heatmap_body(
  "log2ER_high", {
  i = which(colnames(highmat[o,o]) == "macrophage")
  x = i/ncol(highmat[o,o])
  grid.lines(c(x, x), c(1, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black", lineend = "square"))
  grid.lines(c(0, x), c(1-x, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black", lineend = "square"))
  #grid.lines(c(0, x), c(1, 1), gp = gpar(lwd = 4, lty = 1, col = "black"))
  #grid.lines(c(0, 0), c(1, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black"))
  }
)
decorate_heatmap_body(
  "log2ER_low", {
    i = which(colnames(lowmat[o,o]) == "macrophage")
    x = i/ncol(lowmat[o,o])
    grid.lines(c(x, x), c(1, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black", lineend = "square"))
    grid.lines(c(0, x), c(1-x, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black", lineend = "square"))
    #grid.lines(c(0, x), c(1, 1), gp = gpar(lwd = 4, lty = 1, col = "black"))
    #grid.lines(c(0, 0), c(1, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black"))
  }
)
#dev.off()

# ** S6b **
#pdf(file = "7-figures/figS6b.pdf", width = 7, height = 5)
draw(dh,heatmap_legend_side = "bottom")
decorate_heatmap_body(
  "log2FC", {
    i = which(colnames(lowmat[o,o]) == "macrophage")
    x = i/ncol(lowmat[o,o])
    grid.lines(c(x, x), c(1, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black", lineend = "square"))
    grid.lines(c(0, x), c(1-x, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black", lineend = "square"))
    #grid.lines(c(0, x), c(1, 1), gp = gpar(lwd = 4, lty = 1, col = "black"))
    #grid.lines(c(0, 0), c(1, 1-x), gp = gpar(lwd = 4, lty = 1, col = "black"))
  }
)
#dev.off()

# E 
tmas$cellgroup <- ifelse(test = tmas$celltype %in% o[1:14], yes = "immune", no = "other")
fig4e_topleft <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "2f24"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = cellgroup)) + 
  scale_color_manual(values = c("immune"="red", "other"="grey")) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA2 FOV24") +
  NoLegend()
fig4e_topleft
#ggsave(filename = "7-figures/fig4e/fig4e_topleft.pdf", device = "pdf", width = 6, height = 6)

fig4e_topsecondfromleft <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "1f42"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = cellgroup)) + 
  scale_color_manual(values = c("immune"="red", "other"="grey")) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA1 FOV42") +
  NoLegend()
fig4e_topsecondfromleft
#ggsave(filename = "7-figures/fig4e/fig4e_topsecondfromleft.pdf", device = "pdf", width = 6, height = 6)

fig4e_topsecondfromright <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "2f16"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = cellgroup)) + 
  scale_color_manual(values = c("immune"="red", "other"="grey")) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA2 FOV16") +
  NoLegend()
fig4e_topsecondfromright
#ggsave(filename = "7-figures/fig4e/fig4e_topsecondfromright.pdf", device = "pdf", width = 6, height = 6)

fig4e_topright <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "1f44"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = cellgroup)) + 
  scale_color_manual(values = c("immune"="red", "other"="grey")) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA1 FOV44") +
  NoLegend()
fig4e_topright
#ggsave(filename = "7-figures/fig4e/fig4e_topright.pdf", device = "pdf", width = 6, height = 6)

fig4e_legend <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "1f44"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = cellgroup)) + 
  scale_color_manual(values = c("immune"="red", "other"="grey")) + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white"), 
        legend.text = element_text(face = "bold", size = 15), legend.title = element_blank(),
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
legend <- cowplot::get_legend(fig4e_legend)
grid.newpage()
grid.draw(legend)
#ggsave(filename = "7-figures/fig4e/fig4e_legend.pdf", device = "pdf", width = 4, height = 12)

# D
immune_types <- c("T CD4 naive", "mDC", "Treg", "mast", "NK", "T CD8 naive", "pDC", "T CD4 memory", "B-cell", "monocyte", "T CD8 memory", "neutrophil", "plasmablast", "macrophage")
immune_coloco_score <- sum(diffmat[immune_types, immune_types])

dimnames(neres$sim)[[1]] <- gsub(pattern = "celltype", replacement = "", x = dimnames(neres$sim)[[1]]) 
dimnames(neres$sim)[[4]] <- 1:1000
sim_diff <- neres$sim[,,"high",]-neres$sim[,,"low",]
simulated_scores <- sim_diff[immune_types,immune_types,] |> apply(MARGIN = 3, FUN = sum)

fig4d <- ggplot(tibble("score" = simulated_scores)) + 
  geom_histogram(mapping = aes(x = score), bins = 20, color = "black", fill = "grey", alpha = 0.5) + 
  geom_vline(xintercept = immune_coloco_score, colour = "red", linetype = "dashed") + 
  annotate(geom = "text", x = 100, y = 130, size = 6, label = paste("p-value = ", mean(simulated_scores < immune_coloco_score), sep = ""), fontface = "bold") + 
  annotate(geom = "text", x = 100, y = 140, size = 6, label = paste("Observed Difference = ", round(immune_coloco_score, digits = 2), sep = ""), fontface = "bold", color = "red") + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 15)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25)) + 
  labs(x = "score difference", y = "empirical incidences", title = "Randomization Testing for Difference\nin Immune Colocalization Score", subtitle = "High versus Low")
fig4d  
#ggsave("7-figures/fig4d.pdf", device = "pdf", width = 8, height = 8)

# ** S6a **
neres <- readRDS("6-spatial-analysis/line1_cell_level_ne_analysis_results_final_testing.RDS")
rownames(neres$obs) <- gsub(pattern = "celltype", replacement = "", rownames(neres$obs))
rownames(neres$pvals) <- gsub(pattern = "celltype", replacement = "", rownames(neres$pvals))
celltype_cols2 <- c("plasmablast"="magenta",
                   "tumor_L1low"="green4", 
                   "tumor_L1high"="lightgreen",
                   "macrophage"="yellow",
                   "endothelial"="red3", 
                   "fibroblast"="grey", 
                   "mast"="gold4",
                   "B-cell"="purple",
                   "neutrophil"="green", 
                   "cholangiocyte"="dodgerblue", 
                   "Treg"="pink3",
                   "mDC"="yellow3",
                   "T CD8 naive"="orange",
                   "T CD4 memory"="brown2",
                   "monocyte"="darkblue",     
                   "T CD8 memory"="beige",
                   "NK"="tan",
                   "pDC"="orange3", 
                   "T CD4 naive"="turquoise", 
                   "erythrocyte"="pink")  
newo <- c("T CD4 naive",   "mDC",           "Treg",          "mast",          "NK",            "T CD8 naive",   "pDC",           "T CD4 memory", 
          "B-cell",        "monocyte",      "T CD8 memory",  "neutrophil",    "plasmablast",   "macrophage",    "erythrocyte",   "endothelial",  
          "fibroblast",    "cholangiocyte", "tumor_L1low", "tumor_L1high")
colnames(neres$obs) <- gsub(pattern = "_", replacement = "_L1", x = colnames(neres$obs))
rownames(neres$obs) <- gsub(pattern = "_", replacement = "_L1", x = rownames(neres$obs))
mat <- neres$obs[newo, c("tumor_L1high", "tumor_L1low")]

rownames(neres$pvals) <- gsub(pattern = "_", replacement = "_L1", rownames(neres$pvals))
ra <- rowAnnotation(pval = anno_text(case_when((neres$pvals[newo,][,"pval"] < 0.05) & (neres$pvals[newo,][,"fdr"] > 0.05) ~ "*", 
                                               (neres$pvals[newo,][,"fdr"] < 0.05) ~ "**", 
                                               T ~ ""), just = 0, which = "row", gp = gpar(fontface = "bold")) 
)
la <- rowAnnotation(celltype = rownames(mat), 
                    col = list("celltype"=celltype_cols2), 
                    gp = gpar(col = "white", lwd = 0.5), show_legend = F, annotation_name_gp = gpar(col = NA))
ta <- HeatmapAnnotation(celltype = c("tumor_L1high", "tumor_L1low"), 
                        col = list("celltype"=c("tumor_L1low"="green4", "tumor_L1high"="lightgreen")), 
                        gp = gpar(col = "white", lwd = 0.5), show_legend = F, annotation_name_gp = gpar(col = NA))

#pdf(file = "7-figures/figS6a.pdf", width = 4, height = 10)
Heatmap(matrix = mat, name = "log2ER", cluster_columns = F, cluster_rows = F, right_annotation = ra, left_annotation = la, bottom_annotation = ta,
        row_names_side = "left", show_column_names = F,
        col = circlize::colorRamp2(breaks = c(-3, -2, -0.25, 0, 0.25, 1), colors = c(scales::muted("blue"),"blue", "lightblue", "white", "red", scales::muted("red"))),
        #col = viridis::viridis(n = 21),
        width = ncol(mat)*unit(5, "mm"),
        height = nrow(mat)*unit(5, "mm"),
        rect_gp = gpar(col = "white", lwd = 0.5), border = T, border_gp = gpar(col = "black"), 
        row_names_gp = gpar(fontface = "bold"), 
        heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"), 
        column_title = "Colocalization with\nLINE1-ORF1 High and Low\nTumor Cells", column_title_gp = gpar(fontface = "bold")
        ) |> 
  draw(heatmap_legend_side = "bottom")
#dev.off()

# B
tumor_cells <- (tmas$celltype == "tumor")
tmas$line1orf1 <- tmas@assays$RNA@data["LINE1-ORF1",]
meta <- tmas@meta.data[tumor_cells,]
meta$patient <- as.character(meta$patient) |> as.factor()
meta <- dplyr::group_by(meta, patient) |> dplyr::mutate(cut = quantile(line1orf1, 2/3))
meta$line1orf1_status <- ifelse(test = (meta$line1orf1 < meta$cut), yes = "L1low", no = "L1high")
meta$celltype <- paste(meta$celltype, meta$line1orf1_status, sep = "_")
tmas$celltype_tmp <- "other"
tmas@meta.data[meta$cell,]$celltype_tmp <- meta$celltype
tmp_cols <- c("tumor_L1high"="firebrick", "tumor_L1low"="dodgerblue", "other"="grey")
fig4b_top <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "1f19"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype_tmp)) + 
  scale_color_manual(values = tmp_cols) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA1 FOV19") +
  NoLegend()
fig4b_top
#ggsave(filename = "7-figures/fig4b_left.pdf", device = "pdf", width = 6, height = 6)

fig4b_legend_top <-  ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "1f19"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype_tmp)) + 
  scale_color_manual(values = tmp_cols) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white"), 
        legend.text = element_text(face = "bold", size = 15), legend.title = element_blank(),
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
fig4b_legend_top
legend <- cowplot::get_legend(fig4b_legend_top)
grid.newpage()
grid.draw(legend)
#ggsave(filename = "7-figures/fig4b_legend.pdf", device = "pdf", width = 4, height = 12)

#tmas@meta.data[tmas$celltype == "mast",]$celltype_tmp <- "mast"
#tmp_cols <- c("tumor_L1high"="firebrick", "tumor_L1low"="dodgerblue", "other"="grey", "mast"="magenta")
fig4b_bottom <- ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "2f72"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype_tmp)) + 
  scale_color_manual(values = tmp_cols) + 
  #scale_size_manual(values = c("tumor_L1high"=1, "tumor_L1low"=1, "other"=1, "mast"=3)) +
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA2 FOV72") +
  NoLegend()
fig4b_bottom
#ggsave(filename = "7-figures/fig4b_right.pdf", device = "pdf", width = 6, height = 6)

fig4b_legend_bottom <-  ggplot() + 
  geom_point(data = tmas@meta.data |> filter(fov == "2f73"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype_tmp)) + 
  scale_color_manual(values = tmp_cols) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white"), 
        legend.text = element_text(face = "bold", size = 15), legend.title = element_blank(),
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
fig4b_legend_bottom
legend <- cowplot::get_legend(fig4b_legend_bottom)
grid.newpage()
grid.draw(legend)
#ggsave(filename = "7-figures/fig4b_legend_bottom.pdf", device = "pdf", width = 4, height = 12)

# A
dimnames(neres$ptobs)[[1]] <- gsub(pattern = "_", replacement = "_L1", dimnames(neres$ptobs)[[1]])
dimnames(neres$ptobs)[[2]] <- gsub(pattern = "_", replacement = "_L1", dimnames(neres$ptobs)[[2]])
dd <- apply(neres$ptobs, 3, identity, simplify = F) |> purrr::map(as.data.frame) 
dd <- mapply(`[<-`, dd, "celltype", value = purrr::map(dd, .f = rownames), SIMPLIFY = F)
dd <- bind_rows(dd, .id = "patient")
dd <- tidyr::pivot_longer(data = dd, cols = c("tumor_L1high", "tumor_L1low"), names_to = "neighbor", values_to = "log2ER")
dd$celltype <- gsub(pattern = "celltype", replacement = "", x = dd$celltype)
dd$fdr <- plyr::mapvalues(x = dd$celltype, from = rownames(neres$pvals), to = neres$pvals$fdr) |> as.numeric()
dd$pval_position <- case_when(dd$celltype == "mast" ~ 1.5, dd$celltype == "tumor_L1high" ~ 3.5, dd$celltype == "tumor_L1low" ~ 3)
ggplot(data = dd |> filter(celltype %in% c("tumor_L1low", "tumor_L1high"))) + 
  geom_line(mapping = aes(x = neighbor, y = log2ER, group = patient), color = "grey") +
  geom_point(mapping = aes(x = neighbor, y = log2ER)) + 
  geom_text(mapping = aes(x = 1.5, y = pval_position+0.25, label = paste("p.adj = ", round(fdr, digits = 4), sep = "")), check_overlap = T, fontface = "bold", size = 5) +
  stat_summary(mapping = aes(x = neighbor, y = log2ER, color = neighbor), geom = "crossbar", fun = mean) + 
  scale_color_manual(values = c("tumor_L1high"="firebrick", "tumor_L1low"="dodgerblue")) +
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(hjust = 1), axis.title.y = element_text(size = 15),
        strip.text = element_text(face = "bold", size = 20)) +
  facet_wrap(.~celltype, scales = "free_y", nrow = 1) + 
  NoLegend()
#ggsave(filename = "7-figures/fig4a.pdf", device = "pdf", width = 6, height = 6)


## FIGURE 5 --------------------------------------------------------------------
ish <- read.csv("3-DE-analysis/ish-patient-data-de-identified.csv", row.names = 1)

# B
ish$group <- ifelse(test = ish$line1_counts_per_square_um > quantile(ish$line1_counts_per_square_um, 2/3), 
                    yes = "LINE1 HIGH", no = "LINE1 LOW")
ptext <- wilcox.test(x = ish[ish$group == "LINE1 HIGH",]$line1_counts_per_square_um, 
                     y = ish[ish$group == "LINE1 LOW",]$line1_counts_per_square_um, 
                     alternative = "greater", conf.level = 0.95)$p.value
ggplot(data = ish) + 
  stat_summary(mapping = aes(x = group, y = line1_counts_per_square_um), geom = "crossbar", fun = "mean") +
  geom_hline(yintercept = quantile(ish$line1_counts_per_square_um, 2/3), linetype = "dashed") +
  geom_jitter(mapping = aes(x = group, y = line1_counts_per_square_um, shape = group, color = group), width = 0.2, height = 0, size = 3) + 
  scale_color_manual(values = c("LINE1 HIGH"="firebrick", "LINE1 LOW"="dodgerblue")) + 
  scale_shape_manual(values = c("LINE1 HIGH"=16, "LINE1 LOW"=15)) + 
  scale_y_continuous(limits = c(0, 0.08), expand = c(0, 0), breaks = seq(0, 0.08, 0.01)) + 
  annotate(geom = "text", label = "p < 0.0001", x = 2, y = 0.06, fontface = "bold", size = 5) +
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  NoLegend() + 
  theme(axis.title.x = element_blank(), plot.title = element_text(size = 20)) + 
  labs(y = "LINE1 RNA Density per µm²", title = "RNA-ISH Groups")
#ggsave(filename = "7-figures/fig5b.pdf", width = 6, height = 5)


# C
ggplot(data = ish) + 
  geom_point(mapping = aes(x = line1_counts_per_square_um, y = HERVH_PerCell), color = "blue", fill = "blue", shape = 24) + 
  geom_smooth(mapping = aes(x = line1_counts_per_square_um, y = HERVH_PerCell), color = "blue", fill = NA, method = "lm", linewidth = 0.5) + 
  stat_smooth(
    mapping = aes(x = line1_counts_per_square_um, y = HERVH_PerCell),
    method = "lm",
    se = T,
    geom = "ribbon",
    fill = NA,          
    color = "black",     
    alpha = 1            
  ) + 
  annotate(geom = "text", label = paste("R² = ", round(cor(ish$line1_counts_per_square_um, ish$HERVH_PerCell)**2, 4), sep = ""), 
           x = 0.06, y = 185, fontface = "bold") +
  annotate(geom = "text", label = "p < 0.0001", # Per David's analysis
           x = 0.06, y = 177, fontface = "bold") +
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(plot.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 0.08), expand = c(0, 0)) + 
  labs(x = "LINE1 RNA Density per µm²", y = "HERVH Counts per Cell", title = "HERVH")
#ggsave(filename = "7-figures/fig5c_left.pdf", width = 6, height = 6)

ggplot(data = ish) + 
  geom_point(mapping = aes(x = line1_counts_per_square_um, y = HERVK_PerCell), color = "red", fill = "red", shape = 22) + 
  geom_smooth(mapping = aes(x = line1_counts_per_square_um, y = HERVK_PerCell), color = "red", fill = NA, method = "lm", linewidth = 0.5) + 
  stat_smooth(
    mapping = aes(x = line1_counts_per_square_um, y = HERVK_PerCell),
    method = "lm",
    se = TRUE,
    geom = "ribbon",
    fill = NA,          
    color = "black",     
    alpha = 1            
  ) + 
  annotate(geom = "text", label = paste("R² = ", round(cor(ish$line1_counts_per_square_um, ish$HERVK_PerCell)**2, 4), sep = ""), 
           x = 0.06, y = 92, fontface = "bold") +
  annotate(geom = "text", label = "p < 0.0001", # Per David's analysis
           x = 0.06, y = 88, fontface = "bold") +
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(plot.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 0.08), expand = c(0, 0)) + 
  labs(x = "LINE1 RNA Density per µm²", y = "HERVK Counts per Cell", title = "HERVK")
#ggsave(filename = "7-figures/fig5c_middle.pdf", width = 6, height = 6)

ggplot(data = ish) + 
  geom_point(mapping = aes(x = line1_counts_per_square_um, y = HSATII_PerCell), color = "green3", fill = "green3", shape = 23) + 
  geom_smooth(mapping = aes(x = line1_counts_per_square_um, y = HSATII_PerCell), color = "green3", fill = NA, method = "lm", linewidth = 0.5) + 
  stat_smooth(
    mapping = aes(x = line1_counts_per_square_um, y = HSATII_PerCell),
    method = "lm",
    se = TRUE,
    geom = "ribbon",
    fill = NA,          
    color = "black",     
    alpha = 1            
  ) + 
  annotate(geom = "text", label = paste("R² = ", round(cor(ish$line1_counts_per_square_um, ish$HSATII_PerCell)**2, 4), sep = ""), 
           x = 0.06, y = 28, fontface = "bold") +
  annotate(geom = "text", label = "p NS", # Per David's analysis
           x = 0.06, y = 26.75, fontface = "bold") +
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(plot.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 0.08), expand = c(0, 0)) + 
  labs(x = "LINE1 RNA Density per µm²", y = "HSATII Counts per Cell", title = "HSATII")
#ggsave(filename = "7-figures/fig5c_right.pdf", width = 6, height = 6)

# D
library(ggsurvfit)
library(survival)
ish$group <- factor(x = ish$group, levels = c("LINE1 LOW", "LINE1 HIGH"))
fit <- survfit2(Surv(five_year_post_surgery_OS, death) ~ group,
                 data = ish) 
ggsurvfit(x = fit, linewidth = 1.5) +
  labs(
    x = "Days from Surgery",
    y = "Probability of Survival"
  ) +
  #add_risktable() + 
  #add_pvalue(location = "caption") +
  annotate(geom = "text", x = 1750, y = 0.5, label = "p = 0.011", fontface = "bold", size = 5) + 
  scale_color_manual(values = c("LINE1 HIGH" = "firebrick", "LINE1 LOW" = "dodgerblue")) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2000)) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.85), legend.text = element_text(face = "bold", size = 15), 
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
#ggsave(filename = "7-figures/fig5d.pdf", width = 8, height = 8)

## SUPPLEMENTARY FIGURES -------------------------------------------------------

# S1 #
tmas@meta.data |> filter(slide == "hcc_tma1") |> 
  ggplot(mapping = aes(x = CenterX_local_px, y = -CenterY_local_px, color = celltype)) + 
  #geom_point(size = 0.025) +
  scattermore::geom_scattermore(pixels = c(250, 250), pointsize = 1.75, interpolate = F) +
  scale_color_manual(values = celltype_cols) + 
  ggthemes::theme_par() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.position = "top", legend.text.position = "top", 
        panel.background = element_rect(fill = "black"), strip.text = element_text(size = 8)) +
  facet_wrap(.~fov, ncol = 10) + 
  coord_fixed() + 
  NoLegend()
#ggsave(filename = "7-figures/figS1.pdf", width = 8.5, height = 11)

DimPlot(tma1, group.by = "celltype", raster = T) + 
  scale_color_manual(values = celltype_cols) +
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(12, units = "pt"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS1_bottom.pdf", width = 6, height = 6)

# S2 #
tmas@meta.data |> filter(slide == "hcc_tma2") |> 
  ggplot(mapping = aes(x = CenterX_local_px, y = -CenterY_local_px, color = celltype)) + 
  #geom_point(size = 0.025) +
  scattermore::geom_scattermore(pixels = c(250, 250), pointsize = 1.75, interpolate = F) +
  scale_color_manual(values = celltype_cols) + 
  ggthemes::theme_par() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.position = "top", legend.text.position = "top", 
        panel.background = element_rect(fill = "black"), strip.text = element_text(size = 8)) +
  facet_wrap(.~fov, ncol = 10) + 
  coord_fixed() + 
  NoLegend()
#ggsave(filename = "7-figures/figS2.pdf", width = 8.5, height = 11)

DimPlot(tma2, group.by = "celltype", raster = T) + 
  scale_color_manual(values = celltype_cols) +
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(12, units = "pt"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS2_bottom.pdf", width = 6, height = 6)

# S3 #
# A 
tmas@meta.data |> filter(slide == "hcc_tma1") |> dplyr::group_by(fov) |> 
  dplyr::summarise(xmin = min(CenterX_global_px), xmax = max(CenterX_global_px), ymin = min(CenterY_global_px), ymax = max(CenterY_global_px)) |> 
  ggplot() + 
  scattermore::geom_scattermore(data = tmas@meta.data |> filter(slide == "hcc_tma1"), mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype), shape = 21, size = 0.1) +
  scale_color_manual(values = celltype_cols) +
  geom_rect(mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = fov), color = "red", fill = NA) +
  ggrepel::geom_text_repel(mapping = aes(x = (xmax+xmin)/2, y = (ymax+ymin)/2, label = fov), color = "white", face = "bold", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf, size = 6) +
  ggthemes::theme_par() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "top", legend.text.position = "top", 
        panel.background = element_rect(fill = "black"), strip.text = element_text(size = 8)) +
  coord_fixed() + 
  NoLegend()
#ggsave(filename = "7-figures/figS3a_right.pdf", width = 6, height = 12, dpi = 200)

DimPlot(tma1, group.by = "patient", raster = T) + 
  scale_color_manual(values = clcols) +
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(10, units = "pt"), legend.position = "inside", legend.position.inside = c(0.8, 0.8),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS3/figS3b_farright.pdf", width = 6, height = 6)

# B 
tmas@meta.data |> filter(slide == "hcc_tma2") |> dplyr::group_by(fov) |> 
  dplyr::summarise(xmin = min(CenterX_global_px), xmax = max(CenterX_global_px), ymin = min(CenterY_global_px), ymax = max(CenterY_global_px)) |> 
  ggplot() + 
  scattermore::geom_scattermore(data = tmas@meta.data |> filter(slide == "hcc_tma2"), mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype), shape = 21, size = 0.1) +
  scale_color_manual(values = celltype_cols) +
  geom_rect(mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = fov), color = "red", fill = NA) +
  ggrepel::geom_text_repel(mapping = aes(x = (xmax+xmin)/2, y = (ymax+ymin)/2, label = fov), color = "white", face = "bold", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf, size = 6) +
  ggthemes::theme_par() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "top", legend.text.position = "top", 
        panel.background = element_rect(fill = "black"), strip.text = element_text(size = 8)) +
  coord_fixed() + 
  NoLegend()
#ggsave(filename = "7-figures/figS3b_right.pdf", width = 6, height = 12, dpi = 250)

DimPlot(tma2, group.by = "patient", raster = T) + 
  scale_color_manual(values = clcols) +
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(10, units = "pt"), legend.position = "inside", legend.position.inside = c(0.15, 0.9),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS3/figS3b_farright.pdf", width = 6, height = 6)

# C
DimPlot(tmas, group.by = "patient", raster = T) + 
  scale_color_manual(values = clcols) +
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(10, units = "pt"), legend.position = "inside", legend.position.inside = c(0.15, 0.8),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS3/figS3c.pdf", width = 6, height = 6)

tumoronly <- tmas |> subset(celltype == "tumor")
tumoronly <- ScaleData(object = tumoronly, assay = "RNA", features = rownames(tumoronly))
tumoronly <- RunPCA(object = tumoronly, assay = "RNA", features = rownames(tumoronly), npcs = 50)
ElbowPlot(tumoronly, ndims = 50)
tumoronly <- RunUMAP(object = tumoronly, dims = 1:20)

clcols <- grDevices::colors()[!grepl(pattern = "(white)|(grey)|(gray)|(light)", x = grDevices::colors())] |> sample(size = 23, replace = F)
names(clcols) <- tumoronly$patient |> unique()
DimPlot(tumoronly, group.by = "patient", raster = T) + 
  scale_color_manual(values = clcols) + 
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(12, units = "pt"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS3/figS3d.pdf", width = 6, height = 6)

# FeaturePlot(tumoronly, features = "LINE1-ORF1", order = T) + 
#   scale_color_viridis_c() + 
#   scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
#   scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
#   coord_fixed() +
#   theme_classic() +
#   theme(
#     axis.text = element_blank(), 
#     axis.title.x = element_text(face = "bold", size = 12),
#     axis.title.y = element_text(face = "bold", size = 12),
#     axis.ticks = element_blank(),
#     axis.title = element_text(hjust = 0.125),
#     axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
#     legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(12, units = "pt"), legend.position = "inside", legend.position.inside = c(0.7, 0.15),
#     plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
#   )

stromaonly <- tmas |> subset(celltype %in% c("fibroblast", "endothelial"))
stromaonly <- ScaleData(object = stromaonly, assay = "RNA", features = rownames(stromaonly))
stromaonly <- RunPCA(object = stromaonly, assay = "RNA", features = rownames(stromaonly), npcs = 50)
ElbowPlot(stromaonly, ndims = 50)
stromaonly <- RunUMAP(object = stromaonly, dims = 1:15)
DimPlot(stromaonly, group.by = "celltype", pt.size = 3, raster = T) + 
  scale_color_manual(values = celltype_cols) +
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(12, units = "pt"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS3/figS3e_bottom.pdf", width = 6, height = 6)

DimPlot(stromaonly, group.by = "patient", pt.size = 3, raster = T) + 
  scale_color_manual(values = clcols) +
  scale_x_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,1]) / 2, guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = max(tumoronly@reductions$umap@cell.embeddings[,2]) / 2, guide = guide_axis(cap = "upper")) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_blank(), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0.125),
    axis.line = element_line(arrow = arrow(length = unit(3, "pt"), type="closed")), 
    legend.text = element_text(face = "bold", size = 12), legend.key.size = unit(12, units = "pt"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
#ggsave(filename = "7-figures/figS3/figS3e_top.pdf", width = 6, height = 6)

# S4 #
ctde <- read.csv("2-celltyping/celltyping_DE_results.csv", row.names = 1)
ctde$tolab <- ifelse((ctde$logFC > 1.5) & (ctde$fdr < 0.05), yes = T, no = F)
ggplot(data = ctde, mapping = aes(x = logFC, y = -log10(fdr), color = tolab)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey")) +
  geom_vline(xintercept = 1.5) + 
  geom_hline(yintercept = -log10(0.05)) +
  geom_point(size = 0.5) + 
  ggrepel::geom_text_repel(data = ctde[ctde$tolab,], 
                           mapping = aes(x = logFC, y = -log10(fdr), label = feature), 
                           color = "red", min.segment.length = 0, box.padding = 0.1, size = 1.5, max.overlaps = 25, segment.size = 0.25) +
  facet_wrap(.~cluster, scales = "free", ncol = 4) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(strip.text = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) + 
  NoLegend()
#ggsave(filename = "7-figures/figS4.pdf", width = 8.5, height = 11)

# S5 #
# A
endores <- read.csv("3-DE-analysis/line1_patient_level_de_analysis_endothelial_results.csv", row.names = 1)
figs5a <- ggplot() + 
  # Lines ######################################################################
  geom_hline(yintercept = -log10(0.05), linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  # NS points ##################################################################
  geom_point(data = endores |> filter(p.adj > 0.05), 
           mapping = aes(x = logFC, y = -log10(p.adj)), color = "grey") +
  # Left points ################################################################ 
  geom_point(data = endores |> filter(p.adj < 0.05 & logFC < 0), 
           mapping = aes(x = logFC, y = -log10(p.adj)), color = "dodgerblue") +
  ggrepel::geom_text_repel(data = filter(endores, (p.adj < 0.05) & (logFC < 0)), 
                           mapping = aes(x = logFC, y = -log10(p.adj), label = target, fontface = "bold"), 
                           size = 3, max.overlaps = 25, segment.size = 0.25, force_pull = 1, min.segment.length = 0, box.padding = 0.85, ylim = c(1, 8), xlim = c(-1.5, -0.05)) +
  # Right points ###############################################################
  geom_point(data = endores |> filter(p.adj < 0.05 & logFC > 0), 
           mapping = aes(x = logFC, y = -log10(p.adj)), color = "firebrick") +
  ggrepel::geom_text_repel(data = filter(endores, (p.adj < 0.05) & (logFC > 0)), 
                           mapping = aes(x = logFC, y = -log10(p.adj), label = target, fontface = "bold"), 
                           size = 3, max.overlaps = 15, segment.size = 0.25, force_pull = 1, min.segment.length = 0, box.padding = 0.85, ylim = c(1, 8), xlim = c(0.05, 1.5)) +
  # Arrows #####################################################################
  annotate(geom = "segment", x = 0, xend = 0.5, y = -0.5, yend = -0.5, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -0.5, label = "Up in high", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -0.5, yend = -0.5, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = -0.6, y = -0.5, label = "Up in low", hjust = 1, fontface = "bold") +
  # Theme ######################################################################
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  scale_y_continuous(limits = c(-1, 8), breaks = 0:8) +
  scale_x_continuous(limits = c(-1.5, 1.5), breaks = c(-1, 0, 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Endothelial Cell DEGs: LINE1-ORF1 High vs Low Patients")
figs5a
#ggsave(filename = "7-figures/figS5a.pdf", width = 8, height = 8)

# D
ish$patient_deid <- factor(x = ish$patient_deid, levels = ish |> arrange(desc(line1_counts_per_square_um)) |> pull(patient_deid))
fig3a <- ggplot() + 
  geom_bar(data = ish, mapping = aes(x = patient_deid, y = line1_counts_per_square_um, fill = group), stat = "identity") + 
  geom_hline(yintercept = quantile(x = ish$line1_counts_per_square_um, probs = 2/3), linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("LINE1 HIGH" = "firebrick", "LINE1 LOW" = "dodgerblue")) + 
  ggthemes::theme_par() +
  ggpubr::labs_pubr() + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.text = element_text(face = "bold"), legend.position = "inside", legend.position.inside = c(0.75, 0.85), legend.title = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8), 
        axis.text.y = element_text(hjust = 1)) + 
  labs(y = "LINE1 RNA Density per µm²", title = "LINE1 High and Low Patients (by RNA-ISH)")
fig3a
#ggsave(filename = "7-figures/figS5d.pdf", device = "pdf", width = 8, height = 6)

# S6 #
# A -- made above
# B -- made above

