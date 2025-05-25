##### LINE1 DE (Patient Level) Methods Exploration #####
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
library(smiDE)
library(lmerTest)
library(lme4)
library(presto)
library(singlecellmethods)
library(DESeq2)

# Plotting
theme_set(theme_bw())
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
               mapping = aes(x = logFC, y = -log10(p.adj)), alpha = 0.75, stroke = 0.1, fill = cols[1], shape = 21, color = "black") +
    geom_point(data = de_results |> filter(p.adj < 0.05 & logFC > 0), 
               mapping = aes(x = logFC, y = -log10(p.adj)), alpha = 0.75, stroke = 0.1, fill = cols[2], shape = 21, color = "black") +
    ggrepel::geom_text_repel(data = filter(de_results, (p.adj < 0.05) & (target %in% targets_to_lbl)), mapping = aes(x = logFC, y = -log10(p.adj), label = target), 
                             size = 2, max.overlaps = 25, segment.size = 0.2, force_pull = 0.5, min.segment.length = 0) +
    geom_text(data = data.frame(x = c(0.5*(min(de_results$logFC)), 0.5*(max(de_results$logFC))), y = c(-0.25, -0.25), l = plot_labs), 
              mapping = aes(x = x, y = y, label = l)) +
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
tmas$patient <- as.character(tmas$patient)
tmas$line1orf1_group <- plyr::mapvalues(x = tmas$patient, from = ptkey$patient_deid, to = ptkey$line1orf1_group)

all(unique(tmas$patient) %in% ptkey$patient_deid)
# [1] TRUE

# Setting up the data
meta <- dplyr::select(tmas@meta.data, cell, core, patient, slide, fov, sdimx, sdimy, celltype, nCount_RNA, nFeature_RNA, line1orf1_group)
meta$patient <- as.character(meta$patient) |> as.factor()
meta$celltype <- as.factor(meta$celltype)
cts <- tmas@assays$RNA@counts
idx <- meta |> dplyr::arrange(patient) |> rownames()
cts <- cts[,idx]
meta <- meta[idx,]

## TUMOR -----------------------------------------------------------------------
# Start with pseudobulk --> DESeq2
tumoridx <- meta[(meta$celltype == "tumor"),] |> pull(cell)
psb <- collapse_counts(counts_mat = cts[,tumoridx], 
                       meta_data = meta[tumoridx,], 
                       varnames = c("patient", "line1orf1_group"), 
                       get_norm = T, 
                       keep_n = 20) # Ensures that each pseudobulk sample is comprised from >=20 cells

mm <- model.matrix(~line1orf1_group, data = psb$meta_data)
high_vs_low <- (mm[psb$meta_data$line1orf1_group == "high",] |> colMeans()) - (mm[psb$meta_data$line1orf1_group == "low",] |> colMeans())
dds <- DESeq2::DESeqDataSetFromMatrix(countData = psb$counts_mat, colData = psb$meta_data, design = mm)
dds <- estimateSizeFactors(object = dds, type = "poscounts", locfunc = genefilter::shorth)
dds_res <- DESeq2::DESeq(object = dds)
plotDispEsts(object = dds_res)
res <- DESeq2::results(dds_res, contrast = high_vs_low)
res %<>% as.data.frame()
res %<>% dplyr::rename(logFC = log2FoldChange, p.adj = padj)
res$target <- rownames(res)
plot_volcano(de_results = res, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"))
colData(dds)$line1orf1_group %<>% as.factor()
plotCounts(dds = dds, gene = "PDCD1", intgroup = "line1orf1_group")
plotCounts(dds = dds, gene = "POU5F1", intgroup = "line1orf1_group")

res <- res |> na.omit()
mat <- counts(object = dds_res, normalized = T)
mat <- mat[res[res$p.adj < 0.05,]$target, psb$meta_data |> arrange(line1orf1_group) |> rownames()] 
mat <- mat |> apply(MARGIN = 1, FUN = scale) |> t()
Heatmap(matrix = mat, cluster_columns = F,  
        name = "Normalized\nExpression\nz-score", 
        width = ncol(mat)*unit(2, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_names_gp = gpar(fontsize = 6), 
        top_annotation = HeatmapAnnotation(group = psb$meta_data |> arrange(line1orf1_group) |> pull(line1orf1_group), 
                                           col = list(group=c("high"="firebrick", "low"="dodgerblue"))), 
        col = colorRamp2(breaks = c(-4, -2, 0, 2, 4), colors = c("blue", "blue", "white", "red", "red"))
        )

# Next, pseudobulk --> limma
lcts <- edgeR::DGEList(counts = psb$counts_mat)
lcts <- edgeR::calcNormFactors(lcts, method = "TMMwsp")
v <- limma::voom(counts = lcts, design = mm, normalize.method = "quantile", plot = T)
lvfit <- limma::lmFit(object = v, design = mm)
lvres <- lvfit |> limma::contrasts.fit(contrasts = high_vs_low) |> limma::eBayes() |> limma::topTable(number = Inf)
lvres %<>% as.data.frame()
lvres %<>% dplyr::rename(p.adj = adj.P.Val)
lvres$target <- rownames(lvres)
plot_volcano(de_results = lvres, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"))

mat <- v$E[lvres[lvres$p.adj < 0.05,]$target, psb$meta_data |> arrange(line1orf1_group) |> rownames()] 
mat <- mat |> apply(MARGIN = 1, FUN = scale) |> t()
Heatmap(matrix = mat, cluster_columns = F, 
        name = "Normalized\nExpression\nz-score", 
        width = ncol(mat)*unit(2, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_names_gp = gpar(fontsize = 6), 
        top_annotation = HeatmapAnnotation(group = psb$meta_data |> arrange(line1orf1_group) |> pull(line1orf1_group), 
                                           col = list(group=c("high"="firebrick", "low"="dodgerblue"))), 
        #col = viridis::viridis(n = 21)
)

# Next, pseudobulk --> edgeR
dgel <- edgeR::DGEList(counts = psb$counts_mat, group = factor(psb$meta_data$line1orf1_group))
dgel <- edgeR::calcNormFactors(object = dgel, method = "TMMwsp")
dgel <- edgeR::estimateDisp(dgel, design = mm, robust = T)
edgeR::plotBCV(dgel)
efit <- edgeR::glmQLFit(y = dgel, design = mm, robust = T)
edgeR::plotQLDisp(efit)
eres <- edgeR::glmQLFTest(glmfit = efit, contrast = c(0, -1))
eres <- edgeR::topTags(eres, n=Inf) |> as.data.frame()
eres <- dplyr::rename(eres, p.adj = FDR)
eres$target <- rownames(eres)
plot_volcano(eres)

# Next, pseudobulk --> glmgampoi
ggpfit <- glmGamPoi::glm_gp(data = psb$counts_mat, design = mm, size_factors = (dgel$samples$lib.size*dgel$samples$norm.factors))
ggpres <- glmGamPoi::test_de(fit = ggpfit, contrast = c(0, -1))
ggpres <- dplyr::rename(ggpres, p.adj = adj_pval, logFC = lfc, target = name)
plot_volcano(ggpres)

# Next, nebula on the cell level (not pseudobulk) data 
# It is important to have accurate normalization factors in the offset term here.
# DESeq2 and edgeR can calculate these, so I will try both.
# Also, modeling on all cells can allow us to estimate variances more effectively for each gene.
mm <- model.matrix(~celltype+line1orf1_group+celltype:line1orf1_group, data = meta)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = meta, design = mm)
dds <- DESeq2::estimateSizeFactors(object = dds, type = "poscounts")
eff <- (dds@colData$sizeFactor*colSums(cts))

nfs <- edgeR::calcNormFactors(object = cts, method = "TMMwsp")
effs <- (nfs*colSums(cts))

plot(eff, effs)
plot(colSums(cts), eff)
plot(colSums(cts), effs)

# nfit <- nebula::nebula(count = cts, id = meta$patient, pred = mm, covariance = T, offset = colSums(cts), ncore = 10)
# saveRDS(nfit, "3-DE-analysis/nebula_model_patient_level_line1_de_analysis.RDS")
# nfit <- nebula::nebula(count = cts, id = meta$patient, pred = mm, covariance = T, offset = eff, ncore = 10)
# saveRDS(nfit, "3-DE-analysis/nebula_model_patient_level_line1_de_analysis_poscounts.RDS")
# nfit <- nebula::nebula(count = cts, id = meta$patient, pred = mm, covariance = T, offset = effs, ncore = 10)
# saveRDS(nfit, "3-DE-analysis/nebula_model_patient_level_line1_de_analysis_TMMwsp.RDS")

nfit <- readRDS("3-DE-analysis/nebula_model_patient_level_line1_de_analysis.RDS")
nfit_poscounts <- readRDS("3-DE-analysis/nebula_model_patient_level_line1_de_analysis_poscounts.RDS")
nfit_tmmwsp <- readRDS("3-DE-analysis/nebula_model_patient_level_line1_de_analysis_TMMwsp.RDS")

tumor_high_vs_low <- mm[(meta$celltype == "tumor") & (meta$line1orf1_group == "high"), ] |> colMeans() - 
  mm[(meta$celltype == "tumor") & (meta$line1orf1_group == "low"), ] |> colMeans()

tumorout <- local_wald_test(nfit = nfit, .contr = tumor_high_vs_low)
tumorout$convergence <- nfit$convergence 
tumorout <- tumorout[tumorout$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
tumorout$p.adj <- p.adjust(tumorout$p.value, method = "BH")

p1 <- plot_volcano(de_results = tumorout, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"), plot_title = "Library Size")

tumorout_poscounts <- local_wald_test(nfit = nfit_poscounts, .contr = tumor_high_vs_low)
tumorout_poscounts$convergence <- nfit$convergence 
tumorout_poscounts <- tumorout_poscounts[tumorout_poscounts$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
tumorout_poscounts$p.adj <- p.adjust(tumorout_poscounts$p.value, method = "BH")

p2 <- plot_volcano(de_results = tumorout_poscounts, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"), plot_title = "DESeq2: poscounts")

tumorout_tmmwsp <- local_wald_test(nfit = nfit_tmmwsp, .contr = tumor_high_vs_low)
tumorout_tmmwsp$convergence <- nfit$convergence 
tumorout_tmmwsp <- tumorout_tmmwsp[tumorout_tmmwsp$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
tumorout_tmmwsp$p.adj <- p.adjust(tumorout_tmmwsp$p.value, method = "BH")

p3 <- plot_volcano(de_results = tumorout_tmmwsp, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"), plot_title = "edgeR: TMMwsp")

cowplot::plot_grid(plotlist = list(p1, p2, p3), nrow = 1, labels = "AUTO")

mat <- psb$exprs_norm
mat <- mat[tumorout_poscounts |> filter(p.adj < 0.05) |> pull(target), psb$meta_data |> arrange(line1orf1_group) |> rownames()] 
mat <- mat |> apply(MARGIN = 1, FUN = scale) |> t()
Heatmap(matrix = mat, cluster_columns = F, 
        name = "Normalized\nExpression\nz-score", 
        width = ncol(mat)*unit(2, "mm"),
        height = nrow(mat)*unit(0.9, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_names_gp = gpar(fontsize = 3), 
        top_annotation = HeatmapAnnotation(group = psb$meta_data |> arrange(line1orf1_group) |> pull(line1orf1_group), 
                                           col = list(group=c("high"="firebrick", "low"="dodgerblue"))), 
        col = colorRamp2(breaks = c(-4, -2, 0, 2, 4), colors = c("blue", "blue", "white", "red", "red"))
)

# Checking the size factors
tibble(nfs = dds@colData$sizeFactor*colSums(cts), group = dds@colData$line1orf1_group) |> ggplot() + 
  geom_histogram(mapping = aes(x = nfs, fill = group), bins = 25)
tibble(nfs = dds@colData$sizeFactor*colSums(cts), group = dds@colData$line1orf1_group) |> ggplot() + 
  geom_boxplot(mapping = aes(y = nfs, fill = group))
tibble(nfs = nfs, group = dds@colData$line1orf1_group) |> ggplot() + 
  geom_histogram(mapping = aes(x = nfs, fill = group), bins = 25)
tibble(nfs = nfs, group = dds@colData$line1orf1_group) |> ggplot() + 
  geom_boxplot(mapping = aes(y = nfs, fill = group))
tibble(nfs = colSums(cts), group = dds@colData$line1orf1_group) |> ggplot() + 
  geom_histogram(mapping = aes(x = nfs, fill = group), bins = 25)
tibble(nfs = colSums(cts), group = dds@colData$line1orf1_group) |> ggplot() + 
  geom_boxplot(mapping = aes(y = nfs, fill = group))

# Looking at some other cell types
macro_high_vs_low <- mm[(meta$celltype == "macrophage") & (meta$line1orf1_group == "high"), ] |> colMeans() - 
  mm[(meta$celltype == "macrophage") & (meta$line1orf1_group == "low"), ] |> colMeans()

macroout <- local_wald_test(nfit = nfit_poscounts, .contr = macro_high_vs_low)
macroout$convergence <- nfit$convergence 
macroout <- macroout[macroout$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
macroout$p.adj <- p.adjust(macroout$p.value, method = "BH")

plot_volcano(de_results = macroout, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"))

endo_high_vs_low <- mm[(meta$celltype == "endothelial") & (meta$line1orf1_group == "high"), ] |> colMeans() - 
  mm[(meta$celltype == "endothelial") & (meta$line1orf1_group == "low"), ] |> colMeans()

endoout <- local_wald_test(nfit = nfit_poscounts, .contr = endo_high_vs_low)
endoout$convergence <- nfit$convergence 
endoout <- endoout[endoout$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
endoout$p.adj <- p.adjust(endoout$p.value, method = "BH")

plot_volcano(de_results = endoout, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"), plot_title = "Endothelial") +
  ggthemes::theme_par() +
  ggpubr::labs_pubr()

library(ggprism)
plot_volcano(de_results = endoout, plot_labs = c("low", "high"), cols = c("dodgerblue", "firebrick"), plot_title = "Endothelial") +
  ggprism::theme_prism(base_size = 10) + 
  guides(y = "prism_offset", x = "prism_offset")

# Trying out some GSEA
library(msigdbr)
hallmarkgenesets <- msigdbr(species = "Homo sapiens", category = "H") |> 
  dplyr::select(gs_name, gene_symbol)
c5genesets <- msigdbr(species = "Homo sapiens", category = "C5") |> 
  dplyr::select(gs_name, gene_symbol)
c8genesets <- msigdbr(species = "Homo sapiens", category = "C8") |> 
  dplyr::select(gs_name, gene_symbol)

genelist <- lvres |> arrange(desc(t)) |> pull(t)
names(genelist) <- lvres |> arrange(desc(t)) |> pull(target)
tum_hall <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = hallmarkgenesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_hall, "NES")
enrichplot::gseaplot(tum_hall, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

tum_c5 <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = c5genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_c5, "NES", showCategory = 40)
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_NEUROENDOCRINE_CELL_DIFFERENTIATION")
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_WOUND_HEALING")
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_RESPONSE_TO_WOUNDING")

genelist <- res |> arrange(desc(stat)) |> pull(stat)
names(genelist) <- res |> arrange(desc(stat)) |> pull(target)
tum_hall <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = hallmarkgenesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_hall, "NES")
enrichplot::gseaplot(tum_hall, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

tum_c5 <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = c5genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_c5, "NES", showCategory = 40)
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_NEUROENDOCRINE_CELL_DIFFERENTIATION")
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_WOUND_HEALING")
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_RESPONSE_TO_WOUNDING")

tum_c8 <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = c8genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_c8, "NES", showCategory = 10)

genelist <- tumorout_poscounts |> arrange(desc(logFC)) |> pull(logFC)
names(genelist) <- tumorout_poscounts |> arrange(desc(logFC)) |> pull(target)
tum_hall <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = hallmarkgenesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_hall, "NES")
enrichplot::gseaplot(tum_hall, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

genelist <- tumorout_tmmwsp |> arrange(desc(logFC)) |> pull(logFC)
names(genelist) <- tumorout_tmmwsp |> arrange(desc(logFC)) |> pull(target)
tum_hall <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = hallmarkgenesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_hall, "NES")
enrichplot::gseaplot(tum_hall, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

tum_c5 <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = c5genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_c5, "NES", showCategory = 40)
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_NEUROENDOCRINE_CELL_DIFFERENTIATION")
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_WOUND_HEALING")
enrichplot::gseaplot(tum_c5, geneSetID = "GOBP_RESPONSE_TO_WOUNDING")

genelist <- endoout |> arrange(desc(logFC)) |> pull(logFC)
names(genelist) <- endoout |> arrange(desc(logFC)) |> pull(target)
endo_hall <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = hallmarkgenesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(tum_hall, "NES")
enrichplot::gseaplot(tum_hall, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

endo_c5 <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = c5genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(endo_c5, "NES", showCategory = 40)

genelist <- macroout |> arrange(desc(logFC)) |> pull(logFC)
names(genelist) <- macroout |> arrange(desc(logFC)) |> pull(target)
macro_c5 <- clusterProfiler::GSEA(gene = genelist, TERM2GENE = c5genesets, pvalueCutoff = 0.05, pAdjustMethod = "BH")
enrichplot::dotplot(macro_c5, "NES", showCategory = 40)

# Decision: go with the NEBULA result that uses poscounts size factors. 
# Vizualize with volcano plot and dot plot. There will be a dot for each patient, colored by average expression.
# Flag genes that are possibly contaminated. 
# Use theme_par with ggpubr labels. 
# Re-fit model only for the cell types with >1000 cells. 



