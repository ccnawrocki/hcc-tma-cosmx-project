##### LINE1 DE Analysis (Patient Level) #####
## Cole Nawrocki ##

# Environment
rm(list = ls())
.libPaths()

# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project-final/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"    

# Packages
library(magrittr)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)
library(nebula)
library(presto)
library(singlecellmethods)
library(DESeq2)

# Plotting
theme_set(ggthemes::theme_par())
theme_update( 
  plot.title = element_text(hjust = 0.5), 
  panel.grid.minor.y = element_blank(), 
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)

# Functions
local_wald_test <- function(nfit, .contr) {
  
  # Get lfc values
  lfcs <- rowSums(sweep(x = nfit$summary[, 1:length(.contr)], MARGIN = 2, STATS = .contr, FUN = "*"))
  
  # Get p-values
  pvals <- array(data = NA, dim = c(length(lfcs)))
  for (i in 1:length(pvals)) {
    cov <- matrix(NA, length(.contr), length(.contr))
    cov[lower.tri(cov, diag=T)] <- as.numeric(nfit$covariance[i,])
    cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
    pvals[i] <- pchisq(lfcs[i]^2/(t(.contr)%*%cov%*%(.contr)),1,lower.tail=FALSE)
  }
  # Construct output
  out <- data.frame("logFC" = lfcs, 
                    "p.value" = pvals, 
                    "target" = nfit$summary$gene)
  
  return(out)
  
}

plot_volcano <- function(de_results, plot_title = NULL, plot_labs = c(NA, NA), targets_to_lbl = NULL, cols = c("dodgerblue", "red3")) {
  
  if (is.null(targets_to_lbl)) {
    targets_to_lbl <- de_results$target
  }
  
  vp <- ggplot() + 
    geom_hline(yintercept = -log10(0.05), linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    geom_point(data = de_results |> filter(p.adj > 0.05), 
               mapping = aes(x = logFC, y = -log10(p.adj)), alpha = 0.75, stroke = 0.1, color = "grey") +
    geom_point(data = de_results |> filter(p.adj < 0.05 & logFC < 0), 
               mapping = aes(x = logFC, y = -log10(p.adj)), alpha = 0.75, stroke = 0.1, fill = cols[1], color = "black", shape = 21) +
    geom_point(data = de_results |> filter(p.adj < 0.05 & logFC > 0), 
               mapping = aes(x = logFC, y = -log10(p.adj)), alpha = 0.75, stroke = 0.1, fill = cols[2], color = "black", shape = 21) +
    ggrepel::geom_text_repel(data = filter(de_results, (p.adj < 0.05) & (target %in% targets_to_lbl)), mapping = aes(x = logFC, y = -log10(p.adj), label = target), 
                             size = 3, max.overlaps = 15, segment.size = 0.5, force_pull = 0.5, min.segment.length = 0) +
    geom_text(data = data.frame(x = c(0.5*(min(de_results$logFC)), 0.5*(max(de_results$logFC))), y = c(-0.25, -0.25), l = plot_labs), 
              mapping = aes(x = x, y = y, label = l), fontface = "bold", size = 3) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
    labs(title = plot_title) + 
    guides(shape = guide_legend(override.aes = list(size = 5, fill = "black")))
  
  return(vp)
  
}

# Data
tmas <- LoadSeuratRds("hcc_tmas_final.RDS")
ptkey <- read.csv("3-DE-analysis/cosmx-patient-data.csv", row.names = 1)

# Adding necessary data for DE analysis
tmas$core <- paste(substr(x = tmas$slide, start = 8, stop = 8), tmas$core, sep = "")
tmas$fov <- paste(substr(x = tmas$slide, start = 8, stop = 8), tmas$fov, sep = "f")
tmas$patient <- as.character(tmas$patient)
tmas$line1orf1_group <- plyr::mapvalues(x = tmas$patient, from = ptkey$patient_deid, to = ptkey$line1orf1_group)

all(unique(tmas$patient) %in% ptkey$patient_deid)
# [1] TRUE

# Flagging contamination for each cell type
meta <- dplyr::select(tmas@meta.data, cell, core, patient, slide, fov, sdimx, sdimy, celltype, nCount_RNA, nFeature_RNA, line1orf1_group)
cts <- tmas@assays$RNA@counts
# neighbors <- InSituCor:::radiusBasedGraph(x = meta$sdimx, y = meta$sdimy, R = 0.05, subset = meta$fov)
# neighborssymm <- 1 * ((neighbors + Matrix::t(neighbors)) != 0)
# rownames(neighborssymm) <- colnames(neighborssymm) <- colnames(cts)
meta <- dplyr::rename(meta, cell_ID = cell)
# contam <- smiDE:::contamination_ratio_metric(assay_matrix = cts,
#                                              metadata = meta,
#                                              adjacency_matrix = neighborssymm,
#                                              cluster_col = "celltype",
#                                              cellid_col = "cell_ID",
#                                              grouping_col = "fov",
#                                              sdimx_col = "sdimx",
#                                              sdimy_col = "sdimy",
#                                              verbose = TRUE)
# contam <- as.data.frame(contam)
# mean_to_neighborhood_ratio <- matrix(NA, length(unique(contam$target)), length(unique(contam$celltype)),
#                                      dimnames = list(unique(contam$target), unique(contam$celltype)))
# for (cell in unique(contam$celltype)){
#   print(cell)
#   mean_to_neighborhood_ratio[, cell] <- (contam$ratio[contam$celltype == cell])[
#     match(rownames(mean_to_neighborhood_ratio), contam$target[contam$celltype == cell])]
# }
# 
# toomuchcontam <- (mean_to_neighborhood_ratio > 2 | is.na(mean_to_neighborhood_ratio))
# saveRDS(toomuchcontam, file = "3-DE-analysis/contamination.RDS")
toomuchcontam <- readRDS("3-DE-analysis/contamination.RDS")
(!toomuchcontam) |> colSums()

# Setting up the data
meta <- meta |> filter(celltype %in% c("tumor", "endothelial", "macrophage")) 
meta$patient <- as.character(meta$patient) |> as.factor()
meta$celltype <- as.factor(meta$celltype)
idx <- meta |> dplyr::arrange(patient) |> rownames()
cts <- cts[,idx]
meta <- meta[idx,]
mm <- model.matrix(~celltype+line1orf1_group+celltype:line1orf1_group, data = meta)

# Estimating size factors: we could try something like the following: 
# dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = meta, design = mm)
# dds <- DESeq2::estimateSizeFactors(object = dds, type = "poscounts", locfunc = genefilter::shorth)
# eff <- (dds@colData$sizeFactor*colSums(dds@assays@data$counts))

# Martin thinks that it is better to stick with the simplest way for now though.
eff <- colSums(cts)

# Fitting the nebula model
# nfit_final <- nebula::nebula(count = cts, id = meta$patient, pred = mm, covariance = T, offset = eff, ncore = 10)
# saveRDS(nfit_final, file = "3-DE-analysis/nebula_model_patient_level_line1_de_analysis_final.RDS")
nfit_final <- readRDS("3-DE-analysis/nebula_model_patient_level_line1_de_analysis_final.RDS")

# Testing for each cell type and visualizing the results
## TUMOR -----------------------------------------------------------------------
tumor_high_vs_low <- mm[(meta$celltype == "tumor") & (meta$line1orf1_group == "high"), ] |> colMeans() - 
  mm[(meta$celltype == "tumor") & (meta$line1orf1_group == "low"), ] |> colMeans()

tumorout <- local_wald_test(nfit = nfit_final, .contr = tumor_high_vs_low)
tumorout$convergence <- nfit_final$convergence 
tumorout <- tumorout[tumorout$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
tumorout$contamination <- plyr::mapvalues(x = tumorout$target, from = rownames(toomuchcontam), to = toomuchcontam[,"tumor"])
tumorout <- filter(tumorout, contamination == F)
tumorout$p.adj <- p.adjust(tumorout$p.value, method = "BH") # FDR
plot_volcano(de_results = tumorout, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"),
             targets_to_lbl = c("LINE1-ORF1", "LINE1-ORF2", "HSATII", "HERVK", "POU5F1", "LEFTY1", "HCAR2/3", "EFNA4", "EFNA5", "MALAT1", "FGG", "GAS6", "GLUL", "APOA1", "APOC1", "PDCD1", "TNF", "TIGIT"),
             plot_title = "Tumor Cells across LINE1-ORF1\nHigh and Low Patients") + 
  ggpubr::labs_pubr()
ggsave(device = "pdf", width = 6, height = 6, filename = "3-DE-analysis/line1_patient_level_de_analysis_tumor_volcano.pdf")

# Getting profiles for visualization
tmp_norm <- tmas@assays$RNA@data[,tmas$celltype == "tumor"]
tmp_meta <- tmas@meta.data |> filter(celltype == "tumor")
tmp_means_mm <- model.matrix(~0+patient, data = tmp_meta)
tmp_means <- (tmp_norm %*% tmp_means_mm)
tmp_means %<>% as.matrix()
tmp_means <- sweep(x = tmp_means, MARGIN = 2, STATS = colSums(tmp_means_mm), FUN = "/")
colnames(tmp_means) <- gsub(pattern = "patient", replacement = "", x = colnames(tmp_means))

# Getting the proportion of cells that are positive for each gene
propexpressing <- tmp_means * NA
for (pt in colnames(propexpressing)) {
  tempmat <- tmp_norm[rownames(propexpressing), tmp_meta |> filter(patient == pt) |> rownames()]
  propexpressing[, pt] <- Matrix::rowMeans(tempmat > 0)
}

tmp_means <- apply(X = tmp_means, MARGIN = 1, FUN = scale) |> t()
colnames(tmp_means) <- colnames(propexpressing)

# Dot plot
df <- data.frame(gene = rep(rownames(tmp_means), ncol(tmp_means)),
                 patient = rep(colnames(tmp_means), each = nrow(tmp_means)),
                 scaled.mean = as.vector(tmp_means),
                 prop.positive = as.vector(propexpressing))
df <- df[is.element(df$gene, c("LINE1-ORF1", "LINE1-ORF2", "HSATII", "HERVK", "POU5F1", "LEFTY1", "HCAR2/3", "EFNA4", "EFNA5", "MALAT1", "FGG", "GAS6", "GLUL", "APOA1", "APOC1", "PDCD1", "TNF", "TIGIT")), ]

ptorder <- ptkey |> arrange(desc(line1orf1_tumor_cpm)) |> pull(patient_deid)
df$patient <- factor(df$patient, levels = ptorder)

ggplot(data = df, mapping = aes(x=patient, y = gene, fill = scaled.mean, size = prop.positive)) + 
  geom_point(shape = 21, color = "black", stroke = 0.5) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = quantile(df$scaled.mean, probs = c(0.01, 0.99)), 
                       oob = scales::squish) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = ptkey, mapping = aes(x = patient_deid, y = n_distinct(df$gene)+1, fill = line1orf1_group, width = 1, height = 1), inherit.aes = F) +
  scale_fill_manual(values = c("low"="dodgerblue", "high"="firebrick")) +
  scale_x_discrete(expand = c(0, 0)) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
ggsave(device = "pdf", width = 10, height = 6, filename = "3-DE-analysis/line1_patient_level_de_analysis_tumor_dotplot.pdf")

# Heatmap
pdf(file = "3-DE-analysis/line1_patient_level_de_analysis_tumor_heatmap.pdf", width = 8, height = 8)
mat <- tmp_means[c("LINE1-ORF1", "LINE1-ORF2", "HSATII", "HERVK", "POU5F1", "LEFTY1", "HCAR2/3", "EFNA4", "EFNA5", "MALAT1", "FGG", "GAS6", "GLUL", "APOA1", "APOC1", "PDCD1", "TNF", "TIGIT"), ptorder]
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
        col = colorRamp2(breaks = c(-4, -2, 0, 2, 4), colors = c("blue", "blue", "white", "red", "red")), 
        heatmap_legend_param = list(direction = "horizontal")
        ) |> 
  draw(heatmap_legend_side = "bottom")
dev.off()

## ENDOTHELIAL CELLS -----------------------------------------------------------
endothelial_high_vs_low <- mm[(meta$celltype == "endothelial") & (meta$line1orf1_group == "high"), ] |> colMeans() - 
  mm[(meta$celltype == "endothelial") & (meta$line1orf1_group == "low"), ] |> colMeans()

endothelialout <- local_wald_test(nfit = nfit_final, .contr = endothelial_high_vs_low)
endothelialout$convergence <- nfit_final$convergence 
endothelialout <- endothelialout[endothelialout$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
endothelialout$contamination <- plyr::mapvalues(x = endothelialout$target, from = rownames(toomuchcontam), to = toomuchcontam[,"endothelial"])
endothelialout <- filter(endothelialout, contamination == F)
endothelialout$p.adj <- p.adjust(endothelialout$p.value, method = "BH") # FDR
plot_volcano(de_results = endothelialout, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"),
             plot_title = "Endothelial Cells across LINE1-ORF1\nHigh and Low Patients") + 
  ggpubr::labs_pubr()
ggsave(device = "pdf", width = 6, height = 6, filename = "3-DE-analysis/line1_patient_level_de_analysis_endothelial_volcano.pdf")

# Getting profiles for visualization
tmp_norm <- tmas@assays$RNA@data[,tmas$celltype == "endothelial"]
tmp_meta <- tmas@meta.data |> filter(celltype == "endothelial")
tmp_means_mm <- model.matrix(~0+patient, data = tmp_meta)
tmp_means <- (tmp_norm %*% tmp_means_mm)
tmp_means %<>% as.matrix()
tmp_means <- sweep(x = tmp_means, MARGIN = 2, STATS = colSums(tmp_means_mm), FUN = "/")
colnames(tmp_means) <- gsub(pattern = "patient", replacement = "", x = colnames(tmp_means))

# Getting the proportion of cells that are positive for each gene
propexpressing <- tmp_means * NA
for (pt in colnames(propexpressing)) {
  tempmat <- tmp_norm[rownames(propexpressing), tmp_meta |> filter(patient == pt) |> rownames()]
  propexpressing[, pt] <- Matrix::rowMeans(tempmat > 0)
}

tmp_means <- apply(X = tmp_means, MARGIN = 1, FUN = scale) |> t()
colnames(tmp_means) <- colnames(propexpressing)

# Dot plot
df <- data.frame(gene = rep(rownames(tmp_means), ncol(tmp_means)),
                 patient = rep(colnames(tmp_means), each = nrow(tmp_means)),
                 scaled.mean = as.vector(tmp_means),
                 prop.positive = as.vector(propexpressing))
df <- df[is.element(df$gene, endothelialout |> filter(p.adj < 0.05) |> pull(target)), ]

ptorder <- ptkey |> arrange(desc(line1orf1_tumor_cpm)) |> pull(patient_deid)
df$patient <- factor(df$patient, levels = ptorder)

ggplot(data = df, mapping = aes(x=patient, y = gene, fill = scaled.mean, size = prop.positive)) + 
  geom_point(shape = 21, color = "black", stroke = 0.5) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = quantile(df$scaled.mean, probs = c(0.01, 0.99)), 
                       oob = scales::squish) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = ptkey, mapping = aes(x = patient_deid, y = n_distinct(df$gene)+1, fill = line1orf1_group, width = 1, height = 1), inherit.aes = F) +
  scale_fill_manual(values = c("low"="dodgerblue", "high"="firebrick")) +
  scale_x_discrete(expand = c(0, 0)) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size = 8))
ggsave(device = "pdf", width = 10, height = 12, filename = "3-DE-analysis/line1_patient_level_de_analysis_endothelial_dotplot.pdf")

# Heatmap
pdf(file = "3-DE-analysis/line1_patient_level_de_analysis_endothelial_heatmap.pdf", width = 8, height = 12)
mat <- tmp_means[endothelialout |> filter(p.adj < 0.05) |> pull(target), ptorder]
Heatmap(matrix = mat, 
        cluster_columns = F, 
        column_split = (ptkey |> arrange(line1orf1_group) |> pull(line1orf1_group)), 
        column_title_gp = gpar(fontface = "bold"),
        name = "Mean Expression", 
        width = ncol(mat)*unit(2.5, "mm"),
        height = nrow(mat)*unit(2.5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
        column_names_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(group = ptkey |> arrange(line1orf1_group) |> pull(line1orf1_group), 
                                           col = list(group = c("low"="dodgerblue", "high"="firebrick")), 
                                           show_legend = F, 
                                           annotation_name_gp = gpar(col = NA)),
        col = colorRamp2(breaks = c(-4, -2, 0, 2, 4), colors = c("blue", "blue", "white", "red", "red")), 
        heatmap_legend_param = list(direction = "horizontal")
) |> 
  draw(heatmap_legend_side = "bottom")
dev.off()

## MACROPHAGES -----------------------------------------------------------------
macrophage_high_vs_low <- mm[(meta$celltype == "macrophage") & (meta$line1orf1_group == "high"), ] |> colMeans() - 
  mm[(meta$celltype == "macrophage") & (meta$line1orf1_group == "low"), ] |> colMeans()

macrophageout <- local_wald_test(nfit = nfit_final, .contr = macrophage_high_vs_low)
macrophageout$convergence <- nfit_final$convergence 
macrophageout <- macrophageout[macrophageout$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
macrophageout$contamination <- plyr::mapvalues(x = macrophageout$target, from = rownames(toomuchcontam), to = toomuchcontam[,"macrophage"])
macrophageout <- filter(macrophageout, contamination == F)
macrophageout$p.adj <- p.adjust(macrophageout$p.value, method = "BH") # FDR
plot_volcano(de_results = macrophageout, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"),
             plot_title = "Macrophages across LINE1-ORF1\nHigh and Low Patients") + 
  ggpubr::labs_pubr()
ggsave(device = "pdf", width = 6, height = 6, filename = "3-DE-analysis/line1_patient_level_de_analysis_macrophage_volcano.pdf")

# Getting profiles for visualization
tmp_norm <- tmas@assays$RNA@data[,tmas$celltype == "macrophage"]
tmp_meta <- tmas@meta.data |> filter(celltype == "macrophage")
tmp_means_mm <- model.matrix(~0+patient, data = tmp_meta)
tmp_means <- (tmp_norm %*% tmp_means_mm)
tmp_means %<>% as.matrix()
tmp_means <- sweep(x = tmp_means, MARGIN = 2, STATS = colSums(tmp_means_mm), FUN = "/")
colnames(tmp_means) <- gsub(pattern = "patient", replacement = "", x = colnames(tmp_means))

# Getting the proportion of cells that are positive for each gene
propexpressing <- tmp_means * NA
for (pt in colnames(propexpressing)) {
  tempmat <- tmp_norm[rownames(propexpressing), tmp_meta |> filter(patient == pt) |> rownames()]
  propexpressing[, pt] <- Matrix::rowMeans(tempmat > 0)
}

tmp_means <- apply(X = tmp_means, MARGIN = 1, FUN = scale) |> t()
colnames(tmp_means) <- colnames(propexpressing)

# Dot plot
df <- data.frame(gene = rep(rownames(tmp_means), ncol(tmp_means)),
                 patient = rep(colnames(tmp_means), each = nrow(tmp_means)),
                 scaled.mean = as.vector(tmp_means),
                 prop.positive = as.vector(propexpressing))
df <- df[is.element(df$gene, macrophageout |> filter(p.adj < 0.05) |> pull(target)), ]

ptorder <- ptkey |> arrange(desc(line1orf1_tumor_cpm)) |> pull(patient_deid)
df$patient <- factor(df$patient, levels = ptorder)

ggplot(data = df, mapping = aes(x=patient, y = gene, fill = scaled.mean, size = prop.positive)) + 
  geom_point(shape = 21, color = "black", stroke = 0.5) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = quantile(df$scaled.mean, probs = c(0.01, 0.99)), 
                       oob = scales::squish) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = ptkey, mapping = aes(x = patient_deid, y = n_distinct(df$gene)+1, fill = line1orf1_group, width = 1, height = 1), inherit.aes = F) +
  scale_fill_manual(values = c("low"="dodgerblue", "high"="firebrick")) +
  scale_x_discrete(expand = c(0, 0)) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size = 8))
ggsave(device = "pdf", width = 10, height = 8, filename = "3-DE-analysis/line1_patient_level_de_analysis_macrophage_dotplot.pdf")

# Heatmap
pdf(file = "3-DE-analysis/line1_patient_level_de_analysis_macrophage_heatmap.pdf", width = 8, height = 8)
mat <- tmp_means[macrophageout |> filter(p.adj < 0.05) |> pull(target), ptorder]
Heatmap(matrix = mat, 
        cluster_columns = F, 
        column_split = (ptkey |> arrange(line1orf1_group) |> pull(line1orf1_group)), 
        column_title_gp = gpar(fontface = "bold"),
        name = "Mean Expression", 
        width = ncol(mat)*unit(2.5, "mm"),
        height = nrow(mat)*unit(2.5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
        column_names_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(group = ptkey |> arrange(line1orf1_group) |> pull(line1orf1_group), 
                                           col = list(group = c("low"="dodgerblue", "high"="firebrick")), 
                                           show_legend = F, 
                                           annotation_name_gp = gpar(col = NA)),
        col = colorRamp2(breaks = c(-4, -2, 0, 2, 4), colors = c("blue", "blue", "white", "red", "red")), 
        heatmap_legend_param = list(direction = "horizontal")
) |> 
  draw(heatmap_legend_side = "bottom")
dev.off()

# Saving all the results
write.csv(x = tumorout, file = "3-DE-analysis/line1_patient_level_de_analysis_tumor_results.csv")
write.csv(x = endothelialout, file = "3-DE-analysis/line1_patient_level_de_analysis_endothelial_results.csv")
write.csv(x = macrophageout, file = "3-DE-analysis/line1_patient_level_de_analysis_macrophage_results.csv")

# ORA for tumor
library(msigdbr)
c8genesets <- msigdbr(species = "Homo sapiens", category = "C8") |> 
  dplyr::select(gs_name, gene_symbol)
c7genesets <- msigdbr(species = "Homo sapiens", category = "C7") |> 
  dplyr::select(gs_name, gene_symbol)
tftgenesets <- msigdbr(species = "Homo sapiens", category = "C3") |> 
  filter(gs_subcat %in% c("TFT:TFT_Legacy", "TFT:GTRD")) |>
  dplyr::select(gs_name, gene_symbol)
hallmarkgenesets <- msigdbr(species = "Homo sapiens", category = "H") |> 
  dplyr::select(gs_name, gene_symbol)

tum_c8 <- clusterProfiler::enricher(gene = tumorout |> filter(p.adj < 0.05 & logFC > 0) |> pull(target), TERM2GENE = c8genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_c8, "Count")

tum_c7 <- clusterProfiler::enricher(gene = tumorout |> filter(p.adj < 0.05 & logFC > 0) |> pull(target), TERM2GENE = c7genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_c7, x = "Count")

tum_tft <- clusterProfiler::enricher(gene = tumorout |> filter(p.adj < 0.05 & logFC > 0) |> pull(target), TERM2GENE = tftgenesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_tft, x = "Count")

tum_hallmark <- clusterProfiler::enricher(gene = tumorout |> filter(p.adj < 0.05 & logFC > 0) |> pull(target), TERM2GENE = hallmarkgenesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_hallmark, x = "Count")

# This one is interesting
tum_bp <- clusterProfiler::enrichGO(gene = tumorout |> filter(p.adj < 0.05 & logFC > 0) |> pull(target), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_bp, x = "Count" ) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(panel.grid.major = element_line(color = "black", linewidth = 0.25))
ggsave(filename = "3-DE-analysis/line1_patient_level_or_analysis_tumor_dotplot_gobp.pdf", device = "pdf", width = 8, height = 8)

# This one is interesting
tum_mf <- clusterProfiler::enrichGO(gene = tumorout |> filter(p.adj < 0.05 & logFC > 0) |> pull(target), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_mf, x = "Count")
enrichplot::dotplot(tum_mf, x = "Count" ) + 
  ggthemes::theme_par() + 
  ggpubr::labs_pubr() + 
  theme(panel.grid.major = element_line(color = "black", linewidth = 0.25))
ggsave(filename = "3-DE-analysis/line1_patient_level_or_analysis_tumor_dotplot_gomf.pdf", device = "pdf", width = 8, height = 8)

## PSEUDO-BULK -----------------------------------------------------------------
# I just want to try this out.

# Setting up the data
tumor_cells <- (tmas$celltype == "tumor")
meta <- tmas@meta.data[tumor_cells,]
meta$patient <- as.character(meta$patient) |> as.factor()
meta$celltype <- as.factor(meta$celltype)
idx <- meta |> dplyr::arrange(patient) |> rownames()
cts <- tmas@assays$RNA@counts[,idx]
meta <- meta[idx,]

# Groups
meta$line1orf1_group <- plyr::mapvalues(x = meta$patient, from = ptkey$patient_deid, to = ptkey$line1orf1_group)

table(meta$patient, meta$line1orf1_group) # Good balance and plenty of cells in each group for pseudobulking

# Pseudo-bulk is likely best for this, as a paired test will perform well. 
tumor_psb <- collapse_counts(counts_mat = cts, meta_data = meta, varnames = c("patient", "line1orf1_group"), get_norm = F)
colSums(tumor_psb$counts_mat) |> hist(breaks = 25)
mean(tumor_psb$counts_mat == 0)

# Using DESeq2, but with upperquartile for finding size factors.
## Reasoning: 
##  - the data is noisy so I want to avoid normalization --> DESeq2 over limma
##  - the total library sizes are very variable, since some patients were much better-sampled on the TMA --> upperquartile
##  - the data does not have many zeros, so neither TMMwsp nor poscounts are ideal --> upperquartile

dds <- DESeqDataSetFromMatrix(countData = tumor_psb$counts_mat, colData = tumor_psb$meta_data, design = ~line1orf1_group)
nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "upperquartile")
sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
sfs <- sfs / exp(mean(log(sfs))) # Scaling by geometric mean, which DESeq2 expects
dds@colData$sizeFactor <- sfs

desres <- DESeq(dds, fitType = "glmGamPoi", test = "LRT", reduced = ~1) # Using glmGamPoi, which is good for low counts
plotDispEsts(desres) # Dispersion estimates are stable-ish

mm <- model.matrix(~line1orf1_group, data = dds@colData)
hvsl <- mm[dds@colData$line1orf1_group == "high",] |> colMeans() - mm[dds@colData$line1orf1_group == "low",] |> colMeans()
destab <- DESeq2::results(object = desres, contrast = hvsl) |> as.data.frame() # Testing our linear contrast with an LRT

destab <- dplyr::rename(destab, logFC = log2FoldChange, p.adj = padj)
destab$target <- rownames(destab)
destab <- destab |> filter(target %in% rownames(toomuchcontam[!toomuchcontam[,"tumor"],])) # Only testing genes that were not super contaminated
destab$p.adj <- p.adjust(p = destab$p.adj, method = "BH")
plot_volcano(destab |> filter(target != "LINE1-ORF1")) # Plotting (omitting LINE1-ORF1)

# Nothing... not exactly sure what to make of this. I think that there is a lot of patient-specific stuff going on. If 
# I do not model on the cell level or with a paired model, everything else gets sort of overshadowed.

