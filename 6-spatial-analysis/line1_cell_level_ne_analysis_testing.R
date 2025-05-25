# Clear and check environment
rm(list = ls())
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project-final/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"

# Packages
library(magrittr) 
library(circlize) 
library(ComplexHeatmap) 
library(presto)
library(singlecellmethods)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

# Save and set seed
# ssseeeddd <- .Random.seed[1:9]
ssseeeddd <- c(10403, 624, 168352879, 1488529380, -1691384427, 879462738, -110616405, 582342640, -544924879)
set.seed(ssseeeddd)

# Default plotting options
theme_set(theme_bw())
theme_update( 
  plot.title = element_text(hjust = 0.5), 
  panel.grid.minor.y = element_blank(), 
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)

# Data
tmas <- LoadSeuratRds("hcc_tmas_final.RDS")

# Patient IDs
tmas$patient <- as.character(tmas$patient)

# LINE1 expression
tmas$line1orf1 <- tmas@assays$RNA@data["LINE1-ORF1",]

# TUMOR ------------------------------------------------------------------------

# Setting up the data
tumor_cells <- (tmas$celltype == "tumor")
meta <- tmas@meta.data[tumor_cells,]
meta$patient <- as.character(meta$patient) |> as.factor()
meta$celltype <- as.factor(meta$celltype)
idx <- meta |> dplyr::arrange(patient) |> rownames()
cts <- tmas@assays$RNA@counts[,idx]
meta <- meta[idx,]

# Groups
meta <- dplyr::group_by(meta, patient) |> dplyr::mutate(cut = quantile(line1orf1, 2/3))
meta$line1orf1_status <- ifelse(test = (meta$line1orf1 < meta$cut), yes = "low", no = "high")

# Looks good
table(meta$patient, meta$line1orf1_status)

# Extending to all data 
META <- tmas@meta.data
META[meta$cell,]$celltype <- meta$line1orf1_status
META[meta$cell,]$celltype <- paste("tumor", META[meta$cell,]$celltype, sep = "_")
META$celltype |> unique()

META |> group_by(celltype) |> tally() |> pull(n) |> sum()
# [1] 158971

# Running NE analysis we used before, but treated high and low tumor cells as different cell types
# Making FOVs globally-defined
META$fov <- paste(substr(x = META$slide, start = 8, stop = 8), META$fov, sep = "f")
META$fov %<>% as.factor()
META$celltype %<>% as.factor()

# Splitting the data by FOV, since FOVs are not contiguous
META <- split(META, as.character(META$fov))

# Using cell coordinates (in mm) to find Delaunay neighbors (these are first tier neighbors)
coords <- lapply(X = META, FUN = select, sdimx , sdimy)
adj <- purrr::map(coords, spatula::getSpatialNeighbors, return_weights = T)
cellids <- lapply(coords, rownames) |> unlist()
adj <- Matrix::bdiag(adj) # Easier to work with as one big matrix
dimnames(adj) <- list(cellids, cellids)
isSymmetric(adj) # This is a bidirectional spatial network
# [1] TRUE

# We will trim any neighbor distances that are 100 um or greater
trimdists <- function(adj_mat) {
  adj_mat@x[adj_mat@x >= 0.1] <- 0
  adj_mat <- Matrix::drop0(adj_mat)
  return(adj_mat)
}
adj <- trimdists(adj)

# We do not need the distances anymore
rmdists <- function(adj_mat) {
  adj_mat@x[adj_mat@x > 0] <- 1
  return(adj_mat)
}
A <- rmdists(adj)

# We can re-combine the metadata from each FOV
META <- bind_rows(META) |> as.data.frame()
META$patient %<>% as.character() %<>% as.factor()
META$celltype %<>% as.character() %<>% as.factor()
rownames(META) <- META$cell

# Need a function to do the work. 
generate_cellcell_profiles <- function(.adj, .cellmeta, .randomize = F, .typestorandomize = NULL) {

  if (.randomize == T) {
    idx <- (.cellmeta$celltype %in% .typestorandomize)
    ctls <- split(.cellmeta[idx,]$celltype, f = as.character(.cellmeta[idx,]$fov))
    random_ctls <- purrr::map2(.x = ctls, .y = lengths(ctls), sample, replace = F) |> unlist()
    .cellmeta[.cellmeta$celltype %in% .typestorandomize,]$celltype <- random_ctls
    .cellmeta$celltype %<>% as.character() %<>% as.factor()
  }
  
  .adj <- .adj[rownames(.cellmeta), rownames(.cellmeta)]
  
  smm <- sparse.model.matrix(~0+celltype, data = .cellmeta)
  cellcts <- (.adj %*% smm)
  rownames(cellcts) <- rownames(smm)
  colnames(cellcts) <- gsub(pattern = "celltype", replacement = "", x = colnames(cellcts))
  cellcts <- cellcts[(rowSums(cellcts) != 0),]
  smm <- smm[rownames(cellcts),]
  psb_cellcts <- crossprod(x = cellcts, y = smm) |> t() |> as.matrix()
  psb_cellprops <- sweep(x = psb_cellcts, MARGIN = 1, STATS = rowSums(psb_cellcts), FUN = "/")
  exp <- rep(x = colSums(psb_cellcts)/sum(.adj), times = ncol(psb_cellcts)) |> matrix(byrow = T, nrow = nrow(psb_cellcts))
  
  exp[which(colSums(exp) == 0),] <- 0 # This corrects for when a cell type does not exist in a sample!
  
  psb_cellprops[is.nan(psb_cellprops)] <- 0
  exp[is.nan(exp)] <- 0
  
  out <- log2((psb_cellprops+0.0001)/(exp+0.0001))
  
  return(out)
}

# Getting the observed co-localization
obs_ne <- generate_cellcell_profiles(.adj = A, .cellmeta = META, .randomize = F)
Heatmap(obs_ne)

# Randomizing which tumor cells are considered high and low 1000 times, recording co-localization each time
sim_ne <- replicate(n = 1000, 
                    expr = generate_cellcell_profiles(.adj = A, .cellmeta = META, .randomize = T, .typestorandomize = c("tumor_high", "tumor_low")), 
                    simplify = F)
sim_ne <- abind::abind(sim_ne, along = 3)

# We only care about tumor high and low for this
obs_ne <- obs_ne[,c("tumor_high", "tumor_low")]
sim_ne <- sim_ne[,c("tumor_high", "tumor_low"),]

# Empirical p-values
obs_diff <- (obs_ne[,"tumor_high"]-obs_ne[,"tumor_low"])
sim_diff <- (sim_ne[,"tumor_high",]-sim_ne[,"tumor_low",])
get_empircial_p <- function(.index) {
  if (obs_diff[.index] < 0) {
    pval <- mean(sim_diff[.index,] < obs_diff[.index])
  }
  if (obs_diff[.index] > 0) {
    pval <- mean(sim_diff[.index,] > obs_diff[.index])
  }
  return(pval)
}
pvals <- matrix(data = NA, nrow = 20, ncol = 1, dimnames = list(names(obs_diff), "pval"))
for (n in rownames(pvals)) {
  pvals[n,] <- get_empircial_p(n)
}
pvals <- cbind(pvals, "bonferroni"=p.adjust(p = pvals[,"pval"]+0.001, method = "bonferroni"))

# Results
ha <- rowAnnotation(pval = anno_text(case_when((pvals[,"pval"] < 0.05) & (pvals[,"bonferroni"] > 0.05) ~ "*", 
                                               (pvals[,"bonferroni"] < 0.05) ~ "**", 
                                               T ~ ""), just = 0.5, which = "row", gp = gpar(fontface = "bold")) 
                    )
Heatmap(matrix = obs_ne, name = "log2ER", cluster_columns = F, cluster_rows = F, left_annotation = ha,
        col = colorRamp2(breaks = c(-3, -1, 0, 1/3, 1), colors = c("black", "grey", "white", "red", scales::muted("red"))),
        width = ncol(obs_ne)*unit(5, "mm"),
        height = nrow(obs_ne)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1))

# Are the high and low tumor cells really auto-localized? Seems like to an extent.
tmp_cols <- c("tumor_high"="green4", "tumor_low"="gold")
ggplot() + 
  geom_point(data = META |> filter(fov == "2f72"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype)) + 
  scale_color_manual(values = tmp_cols) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA2 FOV72") +
  NoLegend()
ggplot() + 
  geom_point(data = META |> filter(fov == "1f19"), 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = celltype)) + 
  scale_color_manual(values = tmp_cols) + 
  theme_void() + 
  coord_fixed() + 
  theme(plot.background = element_rect(fill = "black"), 
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "white")
  ) + 
  labs(title = "TMA1 FOV19") +
  NoLegend()

# Nonetheless, I do not think that this way of testing is appropriate. We should still 
# be blocking by patient. So, I will conduct paired Wilcoxon rank sum tests.

# Doing the same thing for each patient separately. 
ptidxs <- META$cell |> split(f = META$patient)
ptobs <- list()
for (pt in names(ptidxs)) {
  io <- ptidxs[[pt]]
  ptobs[[pt]] <- generate_cellcell_profiles(.adj = A[io,io], .cellmeta = META[io, ], .randomize = F, .typestorandomize = NULL)
}

# Again, we only care about tumor_high and tumor_low as the neighbors for this.
ptobs <- abind::abind(ptobs, along = 3)[,c("tumor_high", "tumor_low"),]

# Performing the tests
highs <- apply(ptobs[,"tumor_high",], 1, identity, simplify = F)
lows <- apply(ptobs[,"tumor_low",], 1, identity, simplify = F)
wtests <- purrr::map2(.x = highs, .y = lows, .f = wilcox.test, paired = T, conf.level = 0.95, exact = F)
pvals <- sapply(X = wtests, FUN = "[[", "p.value")
pvals <- data.frame("pval"=pvals, "fdr"=p.adjust(p = pvals, method = "BH"), "bonferroni"=p.adjust(p = pvals, method = "bonferroni"))

# New results
ha <- rowAnnotation(pval = anno_text(case_when((pvals[,"pval"] < 0.05) & (pvals[,"fdr"] > 0.05) ~ "*", 
                                               (pvals[,"fdr"] < 0.05) ~ "**", 
                                               T ~ ""), just = 0.5, which = "row", gp = gpar(fontface = "bold")) 
)
Heatmap(matrix = obs_ne, name = "log2ER", cluster_columns = F, cluster_rows = F, left_annotation = ha,
        col = colorRamp2(breaks = c(-3, -1, 0, 1/3, 1), colors = c("black", "grey", "white", "red", scales::muted("red"))),
        width = ncol(obs_ne)*unit(5, "mm"),
        height = nrow(obs_ne)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1))

# More plots
dd <- apply(ptobs, 3, identity, simplify = F) |> purrr::map(as.data.frame) 
dd <- mapply(`[<-`, dd, "celltype", value = purrr::map(dd, .f = rownames), SIMPLIFY = F)
dd <- bind_rows(dd, .id = "patient")
dd <- tidyr::pivot_longer(data = dd, cols = c("tumor_high", "tumor_low"), names_to = "neighbor", values_to = "log2ER")
dd$fdr <- plyr::mapvalues(x = dd$celltype, from = rownames(pvals), to = pvals$fdr) |> as.numeric()
ggplot(data = dd) + 
  geom_line(mapping = aes(x = neighbor, y = log2ER, group = patient)) +
  geom_point(mapping = aes(x = neighbor, y = log2ER)) + 
  stat_summary(mapping = aes(x = neighbor, y = log2ER, color = neighbor), geom = "crossbar", fun = median) + 
  stat_summary(mapping = aes(x = 1.5, y = log2ER+0.5*range(log2ER), label = round(fdr, digits = 4)), geom = "text", fun = max) +
  scale_color_manual(values = c("tumor_high"="firebrick", "tumor_low"="dodgerblue")) +
  facet_wrap(.~celltype, scales = "free")

# Saving 
neres <- list("obs" = obs_ne, 
              "sim" = sim_ne, 
              "ptobs" = ptobs, 
              "pvals" = pvals, 
              "notes" = "test = paired wilcoxon rank sum, design = only tumor cells were split into patient-wise high and low line1 groups")
saveRDS(object = neres, file = "6-spatial-analysis/line1_cell_level_ne_analysis_results_final_testing.RDS")


