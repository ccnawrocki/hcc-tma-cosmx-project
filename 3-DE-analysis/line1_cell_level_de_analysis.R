##### LINE1 DE Analysis (Cell Level) #####
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
    ggnewscale::new_scale_fill() +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
    labs(title = plot_title) + 
    guides(shape = guide_legend(override.aes = list(size = 5, fill = "black")))
  
  return(vp)
  
}

# Data
tmas <- LoadSeuratRds("hcc_tmas_final.RDS")
toomuchcontam <- readRDS("3-DE-analysis/contamination.RDS")
(!toomuchcontam) |> colSums()

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
meta <- dplyr::group_by(meta, patient) |> mutate(cut = quantile(line1orf1, 2/3))
meta$line1orf1_status <- ifelse(test = (meta$line1orf1 < meta$cut), yes = "low", no = "high")

table(meta$patient, meta$line1orf1_status) # Good balance and plenty of cells in each group for pseudobulking

# Pseudo-bulk is likely best for this, as a paired test will perform well. 
tumor_psb <- collapse_counts(counts_mat = cts, meta_data = meta, varnames = c("patient", "line1orf1_status"), get_norm = F)
colSums(tumor_psb$counts_mat) |> hist(breaks = 25)
mean(tumor_psb$counts_mat == 0)

# Using DESeq2, but with upperquartile for finding size factors.
## Reasoning: 
##  - the data is noisy so I want to avoid normalization --> DESeq2 over limma
##  - the total library sizes are very variable, since some patients were much better-sampled on the TMA --> upperquartile
##  - the data does not have many zeros, so neither TMMwsp nor poscounts are ideal --> upperquartile

dds <- DESeqDataSetFromMatrix(countData = tumor_psb$counts_mat, colData = tumor_psb$meta_data, design = ~line1orf1_status+patient)
nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "upperquartile")
sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
sfs <- sfs / exp(mean(log(sfs))) # Scaling by geometric mean, which DESeq2 expects
dds@colData$sizeFactor <- sfs

desres <- DESeq(dds, fitType = "glmGamPoi", test = "LRT", reduced = ~patient) # Using glmGamPoi, which is good for low counts
plotDispEsts(desres) # Dispersion estimates are stable

mm <- model.matrix(~line1orf1_status+patient, data = dds@colData)
hvsl <- mm[dds@colData$line1orf1_status == "high",] |> colMeans() - mm[dds@colData$line1orf1_status == "low",] |> colMeans()
destab <- DESeq2::results(object = desres, contrast = hvsl) |> as.data.frame() # Testing our linear contrast with an LRT

destab <- dplyr::rename(destab, logFC = log2FoldChange, p.adj = padj)
destab$target <- rownames(destab)
destab <- destab |> filter(target %in% rownames(toomuchcontam[!toomuchcontam[,"tumor"],])) # Only testing genes that were not super contaminated
destab$p.adj <- p.adjust(p = destab$p.adj, method = "BH")
plot_volcano(destab |> filter(target != "LINE1-ORF1")) # Plotting (omitting LINE1-ORF1)
ggsave("3-DE-analysis/line1_cell_level_de_analysis_tumor_volcano.pdf", width = 6, height = 6)

# Saving results
write.csv(x = destab, file = "3-DE-analysis/line1_cell_level_de_analysis_tumor_results.csv")

# GSEA
library(msigdbr)

# Genesets for liver cell type signatures
c8genesets <- msigdbr(species = "Homo sapiens", category = "C8") |> 
  dplyr::select(gs_name, gene_symbol)
c8genesets %<>% filter(grepl("LIVER", gs_name))

# Canonical pathways
cpgenesets <- msigdbr(species = "Homo sapiens", category = "C2") |> 
  filter(gs_subcat %in% c("CP")) |>
  dplyr::select(gs_name, gene_symbol)

# Hallmark gene sets
hgenesets <- msigdbr(species = "Homo sapiens", category = "H") |> 
  dplyr::select(gs_name, gene_symbol)

# Cancer modules
cmgenesets <- msigdbr(species = "Homo sapiens", category = "C4") |> 
  filter(gs_subcat %in% c("CM")) |>
  dplyr::select(gs_name, gene_symbol)

# BO:BP
bpgenesets <- msigdbr(species = "Homo sapiens", category = "C5") |> 
  filter(gs_subcat %in% c("GO:BP")) |>
  dplyr::select(gs_name, gene_symbol)

# BO:MF
mfgenesets <- msigdbr(species = "Homo sapiens", category = "C5") |> 
  filter(gs_subcat %in% c("GO:MF")) |>
  dplyr::select(gs_name, gene_symbol)

# Defining our gene list (remove the repeats since they should not be in any gene sets)
destab$signed_stat <- destab$stat*sign(destab$logFC)
tumgenelist <- arrange(destab |> filter(!(target %in% c("LINE1-ORF1", "LINE1-ORF2", "HSATII", "HERVK"))), desc(signed_stat)) |> pull(signed_stat)
names(tumgenelist) <- arrange(destab |> filter(!(target %in% c("LINE1-ORF1", "LINE1-ORF2", "HSATII", "HERVK"))), desc(signed_stat)) |> pull(target)

tum_c8 <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = c8genesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_c8@result$GeneRatio <- (substr(x = tum_c8@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
tum_c8@result$direction <- sign(tum_c8@result$NES)
enrichplot::dotplot(tum_c8, x = "NES", showCategory = 10, split="direction")
enrichplot::gseaplot(x = tum_c8, geneSetID = "AIZARANI_LIVER_C14_HEPATOCYTES_2")
enrichplot::gseaplot(x = tum_c8, geneSetID = "DESCARTES_FETAL_LIVER_HEPATOBLASTS")

tum_cp <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = cpgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_cp@result$GeneRatio <- (substr(x = tum_cp@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
tum_cp@result$direction <- sign(tum_cp@result$NES)
enrichplot::dotplot(tum_cp, x = "NES", showCategory = 10, split="direction")
enrichplot::gseaplot(x = tum_cp, geneSetID = "NABA_SECRETED_FACTORS")

tum_h <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = hgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_h@result$GeneRatio <- (substr(x = tum_h@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
tum_h@result$direction <- sign(tum_h@result$NES)
enrichplot::dotplot(tum_h, x = "NES", showCategory = 10, split="direction")
enrichplot::gseaplot(x = tum_h, geneSetID = "HALLMARK_COMPLEMENT")

tum_cm <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = cmgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_cm@result$GeneRatio <- (substr(x = tum_cm@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
tum_cm@result$direction <- sign(tum_cm@result$NES)
enrichplot::dotplot(tum_cm, x = "NES", showCategory = 10, split="direction")
enrichplot::gseaplot(x = tum_cm, geneSetID = "MODULE_114")

tum_bp <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = bpgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_bp@result$GeneRatio <- (substr(x = tum_bp@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
tum_bp@result$direction <- sign(tum_bp@result$NES)
enrichplot::dotplot(tum_bp, x = "NES", showCategory = 10, split="direction")
enrichplot::gseaplot(x = tum_bp, geneSetID = "GOBP_VESICLE_MEDIATED_TRANSPORT")

tum_mf <- clusterProfiler::GSEA(geneList = tumgenelist, TERM2GENE = mfgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_mf@result$GeneRatio <- (substr(x = tum_mf@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
tum_mf@result$direction <- sign(tum_mf@result$NES)
enrichplot::dotplot(tum_mf, x = "NES", showCategory = 10, split="direction")
enrichplot::gseaplot(x = tum_mf, geneSetID = "GOMF_CELL_ADHESION_MOLECULE_BINDING")
enrichplot::gseaplot(x = tum_mf, geneSetID = "GOMF_GROWTH_FACTOR_ACTIVITY")

# I think really the only interesting ones are C8 and GP:MF

# ENDOTHELIAL ------------------------------------------------------------------

# Setting up the data
endo_cells <- (tmas$celltype == "endothelial")
meta <- tmas@meta.data[endo_cells,]
meta$patient <- as.character(meta$patient) |> as.factor()
meta$celltype <- as.factor(meta$celltype)
idx <- meta |> dplyr::arrange(patient) |> rownames()
cts <- tmas@assays$RNA@counts[,idx]
meta <- meta[idx,]

# Groups
meta <- dplyr::group_by(meta, patient) |> mutate(cut = quantile(line1orf1, 2/3))
meta$line1orf1_status <- ifelse(test = (meta$line1orf1 < meta$cut), yes = "low", no = "high")

table(meta$patient, meta$line1orf1_status) # Not great numbers of cells

# Pseudo-bulk is likely best for this, as a paired test will perform well. 
tumor_psb <- collapse_counts(counts_mat = cts, meta_data = meta, varnames = c("patient", "line1orf1_status"), get_norm = F)
colSums(tumor_psb$counts_mat) |> hist(breaks = 25)
mean(tumor_psb$counts_mat == 0)

dds <- DESeqDataSetFromMatrix(countData = tumor_psb$counts_mat, colData = tumor_psb$meta_data, design = ~line1orf1_status+patient)
nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "TMMwsp")
sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
sfs <- sfs / exp(mean(log(sfs))) # Scaling by geometric mean, which DESeq2 expects
dds@colData$sizeFactor <- sfs

desres <- DESeq(dds, fitType = "local")
plotDispEsts(desres) # Dispersion estimates are stable

mm <- model.matrix(~line1orf1_status+patient, data = dds@colData)
hvsl <- mm[dds@colData$line1orf1_status == "high",] |> colMeans() - mm[dds@colData$line1orf1_status == "low",] |> colMeans()
destab <- DESeq2::results(object = desres, contrast = hvsl) |> as.data.frame() # Testing our linear contrast with an LRT

destab <- dplyr::rename(destab, logFC = log2FoldChange, p.adj = padj)
destab$target <- rownames(destab)
destab <- destab |> filter(target %in% rownames(toomuchcontam[!toomuchcontam[,"endothelial"],])) # Only testing genes that were not super contaminated
destab$p.adj <- p.adjust(p = destab$p.adj, method = "BH")
plot_volcano(destab |> filter(target != "LINE1-ORF1")) # Plotting (omitting LINE1-ORF1)
ggsave("3-DE-analysis/line1_cell_level_de_analysis_endothelial_volcano.pdf", width = 6, height = 6)

# Saving results
write.csv(x = destab, file = "3-DE-analysis/line1_cell_level_de_analysis_endothelial_results.csv")

# MACROPHAGE -------------------------------------------------------------------

# Setting up the data
macro_cells <- (tmas$celltype == "macrophage")
meta <- tmas@meta.data[macro_cells,]
meta$patient <- as.character(meta$patient) |> as.factor()
meta$celltype <- as.factor(meta$celltype)
idx <- meta |> dplyr::arrange(patient) |> rownames()
cts <- tmas@assays$RNA@counts[,idx]
meta <- meta[idx,]

# Groups
meta <- dplyr::group_by(meta, patient) |> mutate(cut = quantile(line1orf1, 2/3))
meta$line1orf1_status <- ifelse(test = (meta$line1orf1 < meta$cut), yes = "low", no = "high")

table(meta$patient, meta$line1orf1_status) # Okay balance, but some groups have low cells

# Pseudo-bulk is likely best for this, as a paired test will perform well. 
tumor_psb <- collapse_counts(counts_mat = cts, meta_data = meta, varnames = c("patient", "line1orf1_status"), get_norm = F)
colSums(tumor_psb$counts_mat) |> hist(breaks = 25)
mean(tumor_psb$counts_mat == 0)

dds <- DESeqDataSetFromMatrix(countData = tumor_psb$counts_mat, colData = tumor_psb$meta_data, design = ~line1orf1_status+patient)
nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "TMMwsp")
sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
sfs <- sfs / exp(mean(log(sfs))) # Scaling by geometric mean, which DESeq2 expects
dds@colData$sizeFactor <- sfs

desres <- DESeq(dds, fitType = "local")
plotDispEsts(desres) # Dispersion estimates are stable

mm <- model.matrix(~line1orf1_status+patient, data = dds@colData)
hvsl <- mm[dds@colData$line1orf1_status == "high",] |> colMeans() - mm[dds@colData$line1orf1_status == "low",] |> colMeans()
destab <- DESeq2::results(object = desres, contrast = hvsl) |> as.data.frame() # Testing our linear contrast with a Wald test

destab <- dplyr::rename(destab, logFC = log2FoldChange, p.adj = padj)
destab$target <- rownames(destab)
destab <- destab |> filter(target %in% rownames(toomuchcontam[!toomuchcontam[,"macrophage"],])) # Only testing genes that were not super contaminated
destab$p.adj <- p.adjust(p = destab$p.adj, method = "BH")
plot_volcano(destab |> filter(target != "LINE1-ORF1")) # Plotting (omitting LINE1-ORF1)

ggsave("3-DE-analysis/line1_cell_level_de_analysis_macrophage_volcano.pdf", width = 6, height = 6)

# Saving results
write.csv(x = destab, file = "3-DE-analysis/line1_cell_level_de_analysis_macrophage_results.csv")

# For the cell types with fewer cells, it is challenging to model with pseudobulking. Trying a cell-level mixed model may be better, but the 
# paired model is working well. Sticking with this model is probably best.


