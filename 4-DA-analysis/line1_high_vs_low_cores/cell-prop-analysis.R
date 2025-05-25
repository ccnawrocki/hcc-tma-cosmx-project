## Cell Proportion Analysis ##
# Cole Nawrocki

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
library(patchwork) # Combining plots
library(ggprism) # Themes
library(circlize) # Color ramps
library(ComplexHeatmap) # Heatmaps
library(ggdendro) # Dendrograms
library(ggh4x) # Dendrograms
## Bioinformatics
library(Seurat) # Processing
library(SeuratDisk) # Read/write
library(smiDE) # DE
library(lmerTest) # DE
library(InSituCor) # Spatial calculations
library(presto)
library(singlecellmethods)
library(DESeq2)
library(multcomp)
library(glmnet)

source("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/functions/Spatial_Functions.R")

# Data
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
coremeta <- readRDS("3-DE-analysis/line1_high_vs_low_cores/de_core_meta.RDS")
demeta <- readRDS("3-DE-analysis/line1_high_vs_low_cores/de_meta.RDS")

pdf("line1_differential_abundance_analysis_plots.pdf", width = 8, height = 8)

## All cell types -------------------------------------------------------------- 
cellcts <- demeta |>
  group_by(core, line1orf1_group, celltype, .drop = F) |> 
  mutate(celltype = as.factor(celltype)) |> 
  tally() |> group_by(core, .drop = F) |> mutate(total = sum(n), prop = n/sum(n)) |> 
  as.data.frame()
props_wide <- pivot_wider(data = cellcts, id_cols = celltype, values_from = prop, values_fill = 0, names_from = core) |> 
  as.data.frame()
rownames(props_wide) <- props_wide$celltype
props_wide <- props_wide[,-1]

coretree = hclust(dist(props_wide |> t()))
dend <- coretree |> as.dendrogram()

ddata_x <- dendro_data(dend)
ddata_x$labels$group <- mapvalues(x = ddata_x$labels$label, from = coremeta$core, to = coremeta$line1orf1_group)
coreorder <- coretree$labels[coretree$order]
cellcts$core <- factor(cellcts$core, levels = coreorder)

p1 <- ggplot(cellcts) + 
  geom_bar(mapping = aes(x = core, y = prop, fill = celltype, group = core), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_dendrogram(hclust = coretree, position = "top", labels = NULL)
p1 + geom_text(data=ddata_x$labels,
               aes(label=label, x=x, y=1.075, color = group), angle = 90) + 
  scale_color_manual(values = c("high"="dodgerblue", "low"="lightblue"))

p2 <- ggplot(cellcts) + 
  geom_bar(mapping = aes(x = core, y = prop, fill = celltype, group = core), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) + 
  scale_x_discrete(limits = coremeta$core[order(coremeta$`cpm.LINE1-ORF1`, decreasing = T)]) + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90))
p2 + geom_text(data=coremeta,
               aes(label=line1orf1_group, x=core, y=1.075, color = line1orf1_group), angle = 90) + 
  scale_color_manual(values = c("high"="dodgerblue", "low"="lightblue"))

props <- group_by(cellcts, line1orf1_group) |> mutate(line1total = sum(n))
props_simple <- group_by(props, line1orf1_group, celltype) |> summarise(prop = sum(n)/line1total)
props_simple <- props_simple[!duplicated(props_simple),]

p3 <- ggplot(props_simple) + 
  geom_bar(mapping = aes(x = line1orf1_group, y = prop, fill = celltype), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) +
  ggthemes::theme_few()
p3

# Testing 
cellcts$patient <- mapvalues(x = cellcts$core, from = coremeta$core, to = coremeta$patient)
cellcts$patient %<>% as.factor()
cellcts$line1orf1_group %<>% as.factor()
mm <- model.matrix(~1+line1orf1_group, data = cellcts)
contr <- mm[cellcts$line1orf1_group == "high",] |> colMeans() - mm[cellcts$line1orf1_group == "low",] |> colMeans()

# Modeling the proportions
proptestfunc <- function(ct) {
  testdata <- filter(cellcts, celltype == ct)
  #modout <- glmer(formula = (n/total) ~ 1+line1orf1_group+(1|patient), weights = testdata$total, data = testdata, family = binomial(link = "logit"))
  modout <- glm(formula = (n/total) ~ 1+line1orf1_group, weights = testdata$total, data = testdata, family = binomial(link = "logit"))
  #modout <- lm(formula = rank(prop) ~ 1+line1orf1_group, data = testdata)
  testout <- (multcomp::glht(model = modout, linfct = matrix(contr, nrow = 1, byrow = T),
                             alternative = "two.sided",
                             rhs = 0)) |> summary(test = univariate())
  outs <- data.frame(
    contrast = "high / low",
    celltype = ct,
    pval = testout[["test"]][["pvalues"]][1],
    Estimate = testout[["test"]][["coefficients"]][[1]],
    tstat = testout[["test"]][["tstat"]][[1]],
    sigma = testout[["test"]][["sigma"]][[1]]
  )
  return(outs)
}
res <- map(.x = unique(cellcts$celltype), .f = proptestfunc)

# Results
res <- bind_rows(res)
res$fdr <- p.adjust(p = res$pval, method = "BH")
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# First model (with random effect)
#             contrast     celltype              pval     Estimate        tstat       sigma        fdr   padj
# -------  -----------  --------------  ----------  -----------  -----------  ----------  ---------  -----
# 1...1    high / low   B-cell           0.6542448    0.4115332    0.4478730   0.9188613   0.852963      1
# 1...2    high / low   NK               0.6524628   -0.2226673   -0.4503435   0.4944388   0.852963      1
# 1...3    high / low   T CD4 memory     0.6800132   -0.2597479   -0.4124451   0.6297757   0.852963      1
# 1...4    high / low   T CD8 memory     0.5857325    0.3615037    0.5450306   0.6632723   0.852963      1
# 1...5    high / low   Treg             0.8229689   -0.1236834   -0.2237281   0.5528294   0.852963      1
# 1...6    high / low   cholangiocyte    0.7146297    0.2761073    0.3656454   0.7551231   0.852963      1
# 1...7    high / low   endothelial      0.4121559   -0.2226717   -0.8201056   0.2715159   0.852963      1
# 1...8    high / low   erythrocyte      0.6616314   -0.3534437   -0.4376618   0.8075726   0.852963      1
# 1...9    high / low   fibroblast       0.2178073   -0.7966349   -1.2323797   0.6464200   0.852963      1
# 1...10   high / low   mDC              0.5419931   -0.3633796   -0.6098018   0.5958979   0.852963      1
# 1...11   high / low   macrophage       0.5396231   -0.1938474   -0.6133830   0.3160300   0.852963      1
# 1...12   high / low   mast             0.8529630    0.0974277    0.1853392   0.5256723   0.852963      1
# 1...13   high / low   monocyte         0.5514415   -0.2790012   -0.5956014   0.4684361   0.852963      1
# 1...14   high / low   neutrophil       0.0778255   -0.7838199   -1.7634448   0.4444823   0.852963      1
# 1...15   high / low   pDC              0.8402847   -0.1159258   -0.2015293   0.5752308   0.852963      1
# 1...16   high / low   tumor            0.8233626   -0.1097624   -0.2232222   0.4917182   0.852963      1
# 1...17   high / low   T CD4 naive      0.6506310   -0.3743528   -0.4528858   0.8265944   0.852963      1
# 1...18   high / low   T CD8 naive      0.4628272   -0.5389535   -0.7341994   0.7340696   0.852963      1
# 1...19   high / low   plasmablast      0.4541763    1.1920319    0.7484706   1.5926236   0.852963      1

# Second model (no random effect)
#          contrast     celltype              pval     Estimate        tstat       sigma         fdr        padj
# -------  -----------  --------------  ----------  -----------  -----------  ----------  ----------  ----------
# 1...1    high / low   B-cell           0.0000000    1.9455930    54.217572   0.0358849   0.0000000   0.0000000
# 1...2    high / low   NK               0.0000000   -0.8590771   -11.022783   0.0779365   0.0000000   0.0000000
# 1...3    high / low   T CD4 memory     0.0000000   -0.9616880    -9.442757   0.1018440   0.0000000   0.0000000
# 1...4    high / low   T CD8 memory     0.0000000    0.5743893     8.924837   0.0643585   0.0000000   0.0000000
# 1...5    high / low   Treg             0.0009966   -0.1931211    -3.291474   0.0586731   0.0010520   0.0189361
# 1...6    high / low   cholangiocyte    0.0000000   -1.0519923   -15.619538   0.0673511   0.0000000   0.0000000
# 1...7    high / low   endothelial      0.0000000    0.1884728     8.459321   0.0222799   0.0000000   0.0000000
# 1...8    high / low   erythrocyte      0.0000000   -0.5329408   -13.247781   0.0402287   0.0000000   0.0000000
# 1...9    high / low   fibroblast       0.0000000   -0.6557512   -25.846978   0.0253705   0.0000000   0.0000000
# 1...10   high / low   mDC              0.0000000   -1.0972191   -11.755757   0.0933346   0.0000000   0.0000000
# 1...11   high / low   macrophage       0.0000000   -0.2993954   -11.786768   0.0254010   0.0000000   0.0000000
# 1...12   high / low   mast             0.0000000   -0.3548528    -5.824070   0.0609287   0.0000000   0.0000001
# 1...13   high / low   monocyte         0.0000002   -0.6523205    -5.208478   0.1252421   0.0000003   0.0000036
# 1...14   high / low   neutrophil       0.0001330   -1.0730234    -3.820761   0.2808402   0.0001487   0.0025278
# 1...15   high / low   pDC              0.0000504   -0.4668104    -4.053771   0.1151546   0.0000598   0.0009576
# 1...16   high / low   tumor            0.0000000   -0.2167135   -17.994579   0.0120433   0.0000000   0.0000000
# 1...17   high / low   T CD4 naive      0.8728432   -0.0428615    -0.160048   0.2678039   0.8728432   1.0000000
# 1...18   high / low   T CD8 naive      0.0000005   -1.0760151    -5.016276   0.2145048   0.0000007   0.0000100
# 1...19   high / low   plasmablast      0.0000000    2.4305831    73.985675   0.0328521   0.0000000   0.0000000

# Third model (wilcoxon rank sum test)
# contrast     celltype              pval   Estimate        tstat      sigma         fdr        padj
# -------  -----------  --------------  ----------  ---------  -----------  ---------  ----------  ----------
# 1...1    high / low   B-cell           0.6365968     -1.650   -0.4776646   3.454306   0.7559587   1.0000000
# 1...2    high / low   NK               0.2990842     -3.600   -1.0580312   3.402546   0.6495649   1.0000000
# 1...3    high / low   T CD4 memory     0.3417124     -3.300   -0.9672249   3.411823   0.6495649   1.0000000
# 1...4    high / low   T CD8 memory     0.9316327     -0.300   -0.0865658   3.465571   0.9316327   1.0000000
# 1...5    high / low   Treg             0.3418763     -3.300   -0.9668912   3.413000   0.6495649   1.0000000
# 1...6    high / low   cholangiocyte    0.5912480      1.875    0.5432628   3.451368   0.7489142   1.0000000
# 1...7    high / low   endothelial      0.2240392     -4.200   -1.2433777   3.377896   0.6495649   1.0000000
# 1...8    high / low   erythrocyte      0.5471220     -2.100   -0.6094780   3.445572   0.7425227   1.0000000
# 1...9    high / low   fibroblast       0.0779769     -6.000   -1.8296030   3.279400   0.6495649   1.0000000
# 1...10   high / low   mDC              0.1764130     -4.650   -1.3869071   3.352784   0.6495649   1.0000000
# 1...11   high / low   macrophage       0.4127668     -2.850   -0.8314273   3.427840   0.6535474   1.0000000
# 1...12   high / low   mast             0.7640552     -1.050   -0.3030987   3.464218   0.8539440   1.0000000
# 1...13   high / low   monocyte         0.2789384     -3.750   -1.1041229   3.396361   0.6495649   1.0000000
# 1...14   high / low   neutrophil       0.0495171     -6.600   -2.0530520   3.214726   0.6495649   0.9408252
# 1...15   high / low   pDC              0.4492582     -2.625   -0.7674148   3.420575   0.6566082   1.0000000
# 1...16   high / low   tumor            0.4127668      2.850    0.8314273   3.427840   0.6535474   1.0000000
# 1...17   high / low   T CD4 naive      0.2691172     -3.525   -1.1274602   3.126496   0.6495649   1.0000000
# 1...18   high / low   T CD8 naive      0.2002646     -4.275   -1.3117320   3.259050   0.6495649   1.0000000
# 1...19   high / low   plasmablast      0.8951258      0.450    0.1330240   3.382848   0.9316327   1.0000000

table(coremeta$core, coremeta$patient) |> colSums()
# 1p1 1p10 1p11 1p13 1p15 1p16 1p17  1p3  1p4  1p5  1p6  1p7  1p8  1p9  2p1 2p10  2p3  2p4  2p5  2p6  2p7  2p8  2p9 
# 1    1    1    2    1    1    2    1    1    1    1    1    1    1    3    1    1    1    2    3    1    1    1 

# If we use the mixed model, there is no statistical power because most of the patients only have one core. If we use 
# the more standard model, everything is significant, which is a sign that there is bias in the study design (I think). 
# Looking at the cores on the flow cell image, I think that this is the case. The cores with the immune aggregates in 
# the FOVs are clustering together. TME cannot really be defined on the core level in this study. It should be defined 
# on the FOV level or sub-FOV level.
ggplot(data = res) + 
  geom_bar(mapping = aes(y = celltype, x = Estimate, fill = padj < 0.05), stat = "identity") + 
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

# Proportion of immune cells ---------------------------------------------------
cellcts <- demeta |> 
  filter(celltype %in% c("B-cell", "mDC", "pDC", "macrophage", "neutrophil", "mast", "plasmablast",
                         "Treg", "T CD4 memory", "T CD8 memory", "T CD4 naive", "T CD8 naive")) |> 
  mutate(celltype = as.factor(celltype)) |>
  group_by(core, line1orf1_group, celltype, .drop = F) |> 
  tally() |> group_by(core) |> mutate(total = sum(n), prop = n/sum(n)) |> 
  as.data.frame()
props_wide <- pivot_wider(data = cellcts, id_cols = celltype, values_from = prop, values_fill = 0, names_from = core) |> 
  as.data.frame()
rownames(props_wide) <- props_wide$celltype
props_wide <- props_wide[,-1]

coretree = hclust(dist(props_wide |> t()))
dend <- coretree |> as.dendrogram()

ddata_x <- dendro_data(dend)
ddata_x$labels$group <- mapvalues(x = ddata_x$labels$label, from = coremeta$core, to = coremeta$line1orf1_group)
coreorder <- coretree$labels[coretree$order]
cellcts$core <- factor(cellcts$core, levels = coreorder)

p1 <- ggplot(cellcts) + 
  geom_bar(mapping = aes(x = core, y = prop, fill = celltype, group = core), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_dendrogram(hclust = coretree, position = "top", labels = NULL)
p1 + geom_text(data=ddata_x$labels,
               aes(label=label, x=x, y=1.075, color = group), angle = 90) + 
  scale_color_manual(values = c("high"="dodgerblue", "low"="lightblue"))

p2 <- ggplot(cellcts) + 
  geom_bar(mapping = aes(x = core, y = prop, fill = celltype, group = core), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) + 
  scale_x_discrete(limits = coremeta$core[order(coremeta$`cpm.LINE1-ORF1`, decreasing = T)]) + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90))
p2 + geom_text(data=coremeta,
               aes(label=line1orf1_group, x=core, y=1.075, color = line1orf1_group), angle = 90) + 
  scale_color_manual(values = c("high"="dodgerblue", "low"="lightblue"))

props <- group_by(cellcts, line1orf1_group) |> mutate(line1total = sum(n))
props_simple <- group_by(props, line1orf1_group, celltype) |> summarise(prop = sum(n)/line1total)
props_simple <- props_simple[!duplicated(props_simple),]

p3 <- ggplot(props_simple) + 
  geom_bar(mapping = aes(x = line1orf1_group, y = prop, fill = celltype), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) +
  ggthemes::theme_few()
p3


# Testing 
cellcts$patient <- mapvalues(x = cellcts$core, from = coremeta$core, to = coremeta$patient)
cellcts$patient %<>% as.factor()
cellcts$line1orf1_group %<>% as.factor()
mm <- model.matrix(~line1orf1_group, data = cellcts)
contr <- mm[cellcts$line1orf1_group == "high",] |> colMeans() - mm[cellcts$line1orf1_group == "low",] |> colMeans()

# Proportion modeled, weighted by sample size (total cells in a core)
proptestfunc <- function(ct) {
  testdata <- filter(cellcts, celltype == ct)
  #modout <- glmer(formula = (n/total) ~ 1+line1orf1_group + (1+1|patient), weights = testdata$total, data = testdata, family = binomial(link = "logit"))
  modout <- glm(formula = (n/total) ~ 1+line1orf1_group, weights = testdata$total, data = testdata, family = binomial(link = "logit"))
  testout <- (multcomp::glht(model = modout, linfct = matrix(contr, nrow = 1, byrow = T),
                             alternative = "two.sided",
                             rhs = 0)) |> summary(test = univariate())
  outs <- data.frame(
    contrast = "high / low",
    celltype = ct,
    pval = testout[["test"]][["pvalues"]][1],
    Estimate = testout[["test"]][["coefficients"]][[1]], 
    tstat = testout[["test"]][["tstat"]][[1]], 
    sigma = testout[["test"]][["sigma"]][[1]]
  )
  return(outs)
}

res <- map(.x = c("B-cell", "mDC", "pDC", "macrophage", "neutrophil", "mast", "plasmablast",
                  "Treg", "T CD4 memory", "T CD8 memory", "T CD4 naive", "T CD8 naive"), .f = proptestfunc)
res <- bind_rows(res)
res$fdr <- p.adjust(p = res$pval, method = "BH")
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

ggplot(data = res) + 
  geom_bar(mapping = aes(y = celltype, x = Estimate, fill = pval < 0.05), stat = "identity") + 
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

# Proportion of T cells --------------------------------------------------------
cellcts <- demeta |>
  filter(celltype %in% c("Treg", "T CD4 memory", "T CD8 memory", "T CD4 naive", "T CD8 naive")) |> 
  mutate(celltype = as.factor(celltype)) |> 
  group_by(core, line1orf1_group, celltype, .drop = F) |> 
  tally() |> group_by(core) |> mutate(total = sum(n), prop = n/sum(n))
props_wide <- pivot_wider(data = cellcts, id_cols = celltype, values_from = prop, values_fill = 0, names_from = core) |> 
  as.data.frame()
rownames(props_wide) <- props_wide$celltype
props_wide <- props_wide[,-1]

coretree = hclust(dist(props_wide |> t()))
dend <- coretree |> as.dendrogram()

ddata_x <- dendro_data(dend)
ddata_x$labels$group <- mapvalues(x = ddata_x$labels$label, from = coremeta$core, to = coremeta$line1orf1_group)
coreorder <- coretree$labels[coretree$order]
props$core <- factor(props$core, levels = coreorder)

p1 <- ggplot(cellcts) + 
  geom_bar(mapping = aes(x = core, y = prop, fill = celltype, group = core), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_dendrogram(hclust = coretree, position = "top", labels = NULL)
p1 + geom_text(data=ddata_x$labels,
               aes(label=label, x=x, y=1.075, color = group), angle = 90) + 
  scale_color_manual(values = c("high"="dodgerblue", "low"="lightblue"))

p2 <- ggplot(cellcts) + 
  geom_bar(mapping = aes(x = core, y = prop, fill = celltype, group = core), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) + 
  scale_x_discrete(limits = coremeta$core[order(coremeta$`cpm.LINE1-ORF1`, decreasing = T)]) + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90))
p2 + geom_text(data=coremeta,
               aes(label=line1orf1_group, x=core, y=1.075, color = line1orf1_group), angle = 90) + 
  scale_color_manual(values = c("high"="dodgerblue", "low"="lightblue"))

props <- group_by(cellcts, line1orf1_group) |> mutate(line1total = sum(n))
props_simple <- group_by(props, line1orf1_group, celltype) |> summarise(prop = sum(n)/line1total)
props_simple <- props_simple[!duplicated(props_simple),]

p3 <- ggplot(props_simple) + 
  geom_bar(mapping = aes(x = line1orf1_group, y = prop, fill = celltype), stat = "identity", color = NA) + 
  scale_fill_manual(values = celltype_cols) +
  ggthemes::theme_few()
p3


# Testing 
cellcts$patient <- mapvalues(x = cellcts$core, from = coremeta$core, to = coremeta$patient)
cellcts$patient %<>% as.factor()
cellcts$line1orf1_group %<>% as.factor()
mm <- model.matrix(~1+line1orf1_group, data = cellcts)
contr <- mm[cellcts$line1orf1_group == "high",] |> colMeans() - mm[cellcts$line1orf1_group == "low",] |> colMeans()

# Proportion modeled, weighted by sample size (total cells in a core)
proptestfunc <- function(ct) {
  testdata <- filter(cellcts, celltype == ct)
  #modout <- glmer(formula = (n/total) ~ 1+line1orf1_group + (1+1|patient), weights = testdata$total, data = testdata, family = binomial(link = "logit"))
  modout <- glm(formula = (n/total) ~ 1+line1orf1_group, weights = testdata$total, data = testdata, family = binomial(link = "logit"))
  testout <- (multcomp::glht(model = modout, linfct = matrix(contr, nrow = 1, byrow = T),
                             alternative = "two.sided",
                             rhs = 0)) |> summary(test = univariate())
  outs <- data.frame(
    contrast = "high / low",
    celltype = ct,
    pval = testout[["test"]][["pvalues"]][1],
    Estimate = testout[["test"]][["coefficients"]][[1]], 
    tstat = testout[["test"]][["tstat"]][[1]], 
    sigma = testout[["test"]][["sigma"]][[1]]
  )
  return(outs)
}

res <- map(.x = c("Treg", "T CD4 memory", "T CD8 memory", "T CD4 naive", "T CD8 naive"), .f = proptestfunc)
res <- bind_rows(res)
res$fdr <- p.adjust(p = res$pval, method = "BH")
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

ggplot(data = res) + 
  geom_bar(mapping = aes(y = celltype, x = Estimate, fill = pval <0.05), stat = "identity") + 
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

dev.off()


