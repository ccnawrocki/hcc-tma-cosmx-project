# Clear and check environment
rm(list = ls())
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project-final/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"

# Packages
library(plyr) 
library(magrittr) 
library(tidyverse) 
library(Matrix) 
library(patchwork)
library(ggprism)
library(circlize) 
library(ComplexHeatmap) 
library(presto)
library(singlecellmethods)
library(Seurat)

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
ptkey <- read.csv("3-DE-analysis/cosmx-patient-data-de-identified.csv", row.names = 1)
meta <- LoadSeuratRds("hcc_tmas_final.RDS")@meta.data
meta$line1orf1_group <- plyr::mapvalues(x = meta$patient, from = ptkey$patient_deid, to = ptkey$line1orf1_group)


# Making FOVs globally-defined
meta$fov <- paste(substr(x = meta$slide, start = 8, stop = 8), meta$fov, sep = "f")
meta$fov %<>% as.factor()
meta$celltype %<>% as.factor()

# Splitting the data by FOV, since FOVs are not contiguous
meta <- split(meta, as.character(meta$fov))

# Using cell coordinates (in mm) to find Delaunay neighbors (these are first tier neighbors)
coords <- lapply(X = meta, FUN = select, sdimx , sdimy)
adj <- purrr::map(coords, spatula::getSpatialNeighbors, return_weights = T)
cellids <- lapply(coords, rownames) |> unlist()
adj <- Matrix::bdiag(adj) # Easier to work with as one big matrix
dimnames(adj) <- list(cellids, cellids)
isSymmetric(adj) # This is a bidirectional spatial network
# [1] TRUE

# Neighbor distances
allneighbordists <- adj |> summary() |> pull(x)
quantile(allneighbordists, 0.99)
#         99% 
# 0.05177191 

tibble(dist = allneighbordists) |> ggplot(aes(dist)) + geom_histogram(bins = 100) + geom_vline(xintercept = 0.1, color = "red") # Linear scale
tibble(dist = allneighbordists) |> ggplot(aes(log10(dist))) + geom_histogram(bins = 100) + geom_vline(xintercept = log10(0.1), color = "red") # Log scale

# We will trim any neighbor distances that are 100 um or greater
trimdists <- function(adj_mat) {
  adj_mat@x[adj_mat@x >= 0.1] <- 0
  adj_mat <- Matrix::drop0(adj_mat)
  return(adj_mat)
}
adj <- trimdists(adj)

# Checking that the trimming worked... it did
allneighbordists <- adj |> summary() |> pull(x)
quantile(allneighbordists, 1)
#       100% 
# 0.09990715 

tibble(dist = allneighbordists) |> ggplot(aes(dist)) + geom_histogram(bins = 100) + geom_vline(xintercept = 0.1, color = "red") # Linear scale
tibble(dist = allneighbordists) |> ggplot(aes(log10(dist))) + geom_histogram(bins = 100) + geom_vline(xintercept = log10(0.1), color = "red") # Log scale

# We do not need the distances anymore
rmdists <- function(adj_mat) {
  adj_mat@x[adj_mat@x > 0] <- 1
  return(adj_mat)
}
A <- rmdists(adj)

# Model matrix for cell types
META <- bind_rows(meta)
META$patient %<>% as.character() %<>% as.factor()
smm <- Matrix::sparse.model.matrix(object = ~0+celltype, data = META)

all(rownames(smm) == rownames(A))
# [1] TRUE

# Getting the make-up for each cell's i-niche (neighbor profile)
Avw <- A %*% smm # Matrix math for adding
Avw %<>% as.matrix() %<>% as.data.frame() # Correcting the data structure
rownames(Avw) <- cellids # Putting the cell IDs back
colnames(Avw) <- gsub(pattern = "celltype", replacement = "", x = colnames(Avw))

# There is one cell that has no neighbors, so we won't consider that one
(rowSums(Avw) == 0) |> which()
Avw <- Avw[-((rowSums(Avw) == 0) |> which()),]
META <- META[rownames(Avw),]

# Now, converting to the patient level 
Avw_list <- Avw |> split(f = META$patient)
Avw_list <- map(.x = Avw_list, .f = t)
ptmeta <- META |> split(f = META$patient)
ptsmms <- lapply(X = ptmeta, FUN = Matrix::sparse.model.matrix, object = ~0+celltype)

ptAvw <- purrr::map2(.x = Avw_list, .y = ptsmms, .f = `%*%`)
ptAvw <- map(.x = ptAvw, .f = as.matrix)
ptobs <- map2(.x = ptAvw, .y = map(.x = ptAvw, .f = colSums), .f = sweep, MARGIN = 2, FUN = "/") # Dividing by sum of Kv (total node degree for Ci)

# Above, we have the observed pairwise neighbor profiles for each patient. Some 
# values are NaN, since some patients have 0 instances of a certain cell type. 

# Here, we will construct the analytically-defined expected pairwise neighbor profiles 
# for each patient.
ptA <- list()
for (pt in levels(META$patient)) {
  idx <- META[META$patient == pt,] |> rownames()
  ptA[[pt]] <- A[idx, idx]
  
}
ptE <- map(.x = ptA, .f = sum)
ptkw <- map(.x = ptAvw, .f = rowSums)
ptexp <- map2(.x = ptkw, .y = ptE, .f = `/`) |> map(.f = `/`, 2)
ptexp <- map(.x = ptexp, .f = rep, 19) |> map(.f = matrix, ncol = 19, byrow = F)
ptexp <- purrr::map2(.x = ptexp, .y = lapply(X = ptobs, FUN = rownames), .f = set_rownames)
ptexp <- purrr::map2(.x = ptexp, .y = lapply(X = ptobs, FUN = colnames), .f = set_colnames)

# To deal with NaNs, we will substitute NaN for 0 and use a small pseudo-count when 
# computing enrichment ratio.
transform_for_ER <- function(mat) {
  mat[is.nan(mat)] <- 0
  return(mat+1e-4)
}
ptobs <- map(.x = ptobs, .f = transform_for_ER)
ptexp <- map(.x = ptexp, .f = transform_for_ER)

pters <- map2(.x = ptobs, .y = ptexp, .f = `/`)
log2ers <- map(.x = pters, .f = log2)

# Looking at a couple heatmaps of the results
Heatmap(log2ers$`2p6`)
Heatmap(log2ers$`1p4`)

# I think that empirical estimation of the expected values will be better for this.

# Now, we need to simulate random neighbors 1000 times for each patient to generate
# an expected pairwise neighbor profile for each patient. This will allow us to 
# get pairwise neighborhood enrichment ratios for each patient.

# First, I need to write a function to do the work.
generate_patient_profiles <- function(.adj, .meta, .randomize = F) {
  
  # .adj is a list of adjacency matrices for each contiguous part of the tissue
  # .meta is a list of meta data for each contiguous part of the tissue
  # .randomize is a logical flag that indicates whether the "celltype" column in .meta should be randomized (the .itype values are left alone) 
  
  if (.randomize == T) {
    ct_list <- sapply(X = .meta, FUN = "[[", "celltype")
    rand_ct_list <- purrr::map2(.x = ct_list, .y = lengths(ct_list), .f = sample, replace = F)
    .meta <- mapply(`[<-`, .meta, "celltype", value = rand_ct_list, SIMPLIFY = FALSE)
  }
  
  smms <- lapply(X = .meta, FUN = Matrix::sparse.model.matrix, object = ~0+celltype)
  cellids <- lapply(X = smms, FUN = rownames)
  
  cellcts_list <- purrr::map2(.x = .adj, .y = smms, .f = `%*%`) # Matrix math for adding
  cellcts_list <- purrr::map(.x = cellcts_list, .f = as.matrix) |> purrr::map(.f = as.data.frame) # Correcting the data structure
  cellcts_list <- purrr::map2(.x = cellcts_list, .y = cellids, .f = set_rownames) # Putting the cell IDs back
  cellcts <- bind_rows(cellcts_list) # Making one big data frame
  colnames(cellcts) <- gsub(pattern = "celltype", replacement = "", x = colnames(cellcts))
  cellcts <- cellcts[(rowSums(cellcts) != 0),]
  
  META <- bind_rows(.meta)[rownames(cellcts),]
  META$patient %<>% as.character() %<>% as.factor()
  cellcts_list <- cellcts |> split(f = META$patient)
  cellcts_list <- map(.x = cellcts_list, .f = t)
  ptmeta <- META |> split(f = META$patient)
  ptsmms <- lapply(X = ptmeta, FUN = Matrix::sparse.model.matrix, object = ~0+celltype)
  
  ptcellcts <- purrr::map2(.x = cellcts_list, .y = ptsmms, .f = `%*%`)
  ptcellcts <- map(.x = ptcellcts, .f = as.matrix)
  ptcellprops <- map2(.x = ptcellcts, .y = map(.x = ptcellcts, .f = colSums), .f = sweep, MARGIN = 2, FUN = "/")
  
  
  return(ptcellprops)
  
}

adj <- purrr::map(coords, spatula::getSpatialNeighbors, return_weights = T)
adj <- purrr::map(.x = adj, .f = trimdists)
adj <- purrr::map(.x = adj, .f = rmdists)

# Getting observed and expected via the function
obs <- generate_patient_profiles(.adj = adj, .meta = meta, .randomize = F)
obs <- map(.x = obs, .f = transform_for_ER)
timing <- system.time(expr = {
  expect <- replicate(n = 1000, expr = generate_patient_profiles(.adj = adj, .meta = meta, .randomize = T), simplify = F)
}
)
timing
#    user  system elapsed 
# 233.243  20.523 254.066 

expect <- map(.x = expect, .f = map, transform_for_ER)

# We have 4 dimensions: index, neighbor, patient, and replication
# Right now, we have a list of 1000. Each element is a list of 23. Each of these is a 2d matrix. 
# I want to create a list of 23, with each element being a 3d array with the following dimensions: 
  # neighbor, index, replication

# First create a list of 3d arrays. This is a list of 1000, with each element being a 19x19x23 3d array:
listofcubes <- map(.x = expect, .f = abind::abind, along = 3)

# Convert to a 4d tensor
d4tensor <- abind::abind(listofcubes, along = 4)
dimnames(d4tensor)[[4]] <- 1:1000
dim(d4tensor)
# [1]   19   19   23 1000

# We have all four dimensions accounted for still. For orientation, let's print a corner:
abind::acorn(d4tensor)
# , , 1p1, 1
# 
#               celltypeB-cell celltypecholangiocyte celltypeendothelial celltypeerythrocyte celltypefibroblast
# B-cell                0.0001            0.02867143         0.007792308          0.01548462         0.01128211
# cholangiocyte         0.0501            0.00010000         0.027023077          0.01804872         0.01287955
# endothelial           0.0501            0.10010000         0.046253846          0.04112564         0.05920543
# erythrocyte           0.1501            0.10010000         0.061638462          0.13343333         0.07517987
# fibroblast            0.1751            0.11438571         0.142407692          0.12061282         0.14067508
# macrophage            0.1001            0.11438571         0.169330769          0.18471538         0.21096262

# Making the list of 23 with the intended dimensions.
newlistofcubes <- apply(d4tensor, 3, identity, simplify = F)

# Getting the averages for each patient, which will serve as the expected values
expected <- map(.x = newlistofcubes, .f = apply, MARGIN = c(1,2), FUN = mean)

# Note: we could also have invoked the operations like this for a more stream-lined method:  
# expected <- apply(X = d4tensor, MARGIN = c(1,2,3), FUN = mean) |> apply(MARGIN = 3, FUN = identity, simplify = F)

# Enrichment scores for each patient
ers <- map2(.x = obs, .y = expected, .f = `/`)
ers <- purrr::map(.x = ers, .f = set_colnames, gsub(pattern = "celltype", replacement = "", x = colnames(ers$`1p1`)))

# Order of cell types that makes sense
idx <-  c("Treg", "plasmablast", "B-cell", "T CD8 memory", "T CD8 naive", "T CD4 memory", "T CD4 naive", "NK", "pDC", "mDC", "macrophage")

# All line1-high patients
pdf(file = "6-spatial-analysis/line1-high-patients-NE-split-by-patient.pdf", width = 6, height = 6)
highpts <- ptkey |> filter(line1orf1_group == "high") |> pull(patient_deid)
for (pt in highpts) {
  m <- log2(ers[[pt]]) 
  h <- 
    Heatmap(matrix = m[idx,idx], cluster_columns = F, cluster_rows = F, name = "log2ER", 
            col = circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c(scales::muted("blue"), "white", scales::muted("red"))), 
            width = ncol(m)*unit(5, "mm"),
            height = nrow(m)*unit(5, "mm"), 
            rect_gp = gpar(col = "white", lwd = 0.5), 
            row_title = "Neighbor", column_title = "Index")
  draw(h, column_title = pt)
}
dev.off()

# All line1-low patients
pdf(file = "6-spatial-analysis/line1-low-patients-NE-split-by-patient.pdf", width = 6, height = 6)
lowpts <- ptkey |> filter(line1orf1_group == "low") |> pull(patient_deid)
for (pt in lowpts) {
  m <- log2(ers[[pt]]) 
  h <- 
    Heatmap(matrix = m[idx,idx], cluster_columns = F, cluster_rows = F, name = "log2ER", 
            col = circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c(scales::muted("blue"), "white", scales::muted("red"))), 
            width = ncol(m)*unit(5, "mm"),
            height = nrow(m)*unit(5, "mm"), 
            rect_gp = gpar(col = "white", lwd = 0.5), 
            row_title = "Neighbor", column_title = "Index")
  draw(h, column_title = pt)
}
dev.off()

# Average for high patients and average for low patients: 
highmean <- pters[highpts] |> map(.f = log2) |> abind::abind(along = 3) |> apply(MARGIN = c(1,2), FUN = mean)
colnames(highmean) <- gsub(pattern = "celltype", replacement = "", x = colnames(highmean))
lowmean <- pters[lowpts] |> map(.f = log2) |> abind::abind(along = 3) |> apply(MARGIN = c(1,2), FUN = mean)
colnames(lowmean) <- gsub(pattern = "celltype", replacement = "", x = colnames(lowmean))

idx <-  c("Treg", "plasmablast", "B-cell", "T CD8 memory", "T CD8 naive", "T CD4 memory", "T CD4 naive", "NK", "pDC", "mDC", "macrophage")
o <- hclust(dist(lowmean[idx, idx]))$order
hh <- Heatmap(matrix = highmean[idx,idx][o,o], cluster_columns = F, cluster_rows = F, name = "Mean\nlog2ER", 
        col = circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c(scales::muted("blue"), "white", scales::muted("red"))), 
        width = ncol(m)*unit(5, "mm"),
        height = nrow(m)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_title = "Neighbor", column_title = "Index")
hl <- Heatmap(matrix = lowmean[idx,idx][o,o], cluster_columns = F, cluster_rows = F, name = "Mean\nlog2ER", 
        col = circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c(scales::muted("blue"), "white", scales::muted("red"))), 
        width = ncol(m)*unit(5, "mm"),
        height = nrow(m)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5), 
        row_title = "Neighbor", column_title = "Index")
pdf(file = "6-spatial-analysis/line1-high-and-low-average-NE.pdf", width = 6, height = 6)
draw(hh, column_title = "LINE1-High Patients")
draw(hl, column_title = "LINE1-Low Patients")
dev.off()

# "Immune organization meta-score" will be defined as the sum of the red immune boxes of interest
immune_meta_score <- function(.mat, .idx) {
  .mat[.mat < 0] <- 0
  return(.mat[.idx,.idx] |> sum())
}
scores <- map(.x = ers |> map(.f = log2), .f = immune_meta_score, .idx = idx) 
scores %<>% unlist() %<>% as.data.frame()
scores$group <- plyr::mapvalues(x = rownames(scores), from = ptkey$patient_deid, to = ptkey$line1orf1_group)
colnames(scores) <- c("immune_score", "line1orf1_group")
scores$patient <- rownames(scores)

# Bar plot of the meta-score, split by high and low
ggplot(scores) + 
  geom_boxplot(mapping = aes(x = line1orf1_group, y = immune_score), outliers = F) +
  geom_jitter(mapping = aes(x = line1orf1_group, y = immune_score), height = 0, width = 0.1) + 
  stat_summary(mapping = aes(x = line1orf1_group, y = immune_score, color = line1orf1_group), geom = "crossbar", fun = mean) + 
  scale_color_manual(values = c("high"="firebrick", "low"="dodgerblue"))
ggsave(filename = "6-spatial-analysis/line1_high_and_low_immune_scores_boxplot.pdf", height = 8, width = 6, device = "pdf")

# Wilcoxon rank sum test 
wilcox.test(x = scores |> filter(line1orf1_group == "high") |> pull(immune_score), 
            y = scores |> filter(line1orf1_group == "low") |> pull(immune_score))
# p-value = 0.591

# Bootstrap CIs
mean.boot <- function(d, i) {
  return(mean(d[i]))
}
highboot <- boot::boot(scores |> filter(line1orf1_group == "high") |> pull(immune_score), mean.boot, R = 1000)
boot::boot.ci(highboot, type = "basic")$basic[,c(4,5)]
# 33.93719 81.23806 

lowboot <- boot::boot(scores |> filter(line1orf1_group == "low") |> pull(immune_score), mean.boot, R = 1000)
boot::boot.ci(lowboot, type = "basic")$basic[,c(4,5)]
# 47.50879 91.70717

# These overlap, so we cannot reject the null that the immune score is higher in the low group.

# Question we were trying to answer: is there more immune organization on average in patients with line1-low tumors?
# Answer: Not based on signficance testing. Effect is in that direction though.

# Last idea: inspired from cellcharter 
#   Compute the obs and exp, using the analystical approach, for both high and low. 
#   Mix up the high and low labels 1000 times and compute high and low each time. 
#   Create empirical distribution for high vs low for each cell-cell interaction. 
#   P-value = based on this empirical distribution.
# The key here is that the patients are being mixed up.

# Again, need a function to do the work. 
generate_cellcell_profiles <- function(.adj, .cellmeta, .ptmeta, .factor, .levels, .randomize = F) {
  
  if (.randomize == T) {
    newf <- sample(x = ptkey$line1orf1_group, size = nrow(ptkey), replace = F)
    .ptmeta[[.factor]] <- newf
  }
  
  idx1 <- .ptmeta |> filter(get(.factor) == .levels[1]) |> pull(patient_deid)
  idx2 <- .ptmeta |> filter(get(.factor) == .levels[2]) |> pull(patient_deid)
  m1 <- .cellmeta |> filter(patient %in% idx1)
  m2 <- .cellmeta |> filter(patient %in% idx2)
  A1 <- .adj[rownames(m1), rownames(m1)]
  A2 <- .adj[rownames(m2), rownames(m2)]
  
  smm1 <- sparse.model.matrix(~0+celltype, data = m1)
  cellcts1 <- (A1 %*% smm1)
  rownames(cellcts1) <- rownames(smm1)
  colnames(cellcts1) <- gsub(pattern = "celltype", replacement = "", x = colnames(cellcts1))
  cellcts1 <- cellcts1[(rowSums(cellcts1) != 0),]
  psb_cellcts1 <- crossprod(x = cellcts1, y = smm1) |> t() |> as.matrix()
  psb_cellprops1 <- sweep(x = psb_cellcts1, MARGIN = 1, STATS = rowSums(psb_cellcts1), FUN = "/")
  exp1 <- rep(x = colSums(psb_cellcts1)/sum(A1), times = ncol(psb_cellcts1)) |> matrix(byrow = T, nrow = nrow(psb_cellcts1))
  exp1[which(colSums(exp1) == 0),] <- 0 # This corrects for when a cell type does not exist in a sample!
  
  smm2 <- sparse.model.matrix(~0+celltype, data = m2)
  cellcts2 <- (A2 %*% smm2)
  rownames(cellcts2) <- rownames(smm2)
  colnames(cellcts2) <- gsub(pattern = "celltype", replacement = "", x = colnames(cellcts2))
  cellcts2 <- cellcts2[(rowSums(cellcts2) != 0),]
  psb_cellcts2 <- crossprod(x = cellcts2, y = smm2) |> t() |> as.matrix()
  psb_cellprops2 <- sweep(x = psb_cellcts2, MARGIN = 1, STATS = rowSums(psb_cellcts2), FUN = "/")
  exp2 <- rep(x = colSums(psb_cellcts2)/sum(A2), times = ncol(psb_cellcts2)) |> matrix(byrow = T, nrow = nrow(psb_cellcts2))
  exp2[which(colSums(exp2) == 0),] <- 0 # This corrects for when a cell type does not exist in a sample!
  
  psb_cellprops1[is.nan(psb_cellprops1)] <- 0
  psb_cellprops2[is.nan(psb_cellprops2)] <- 0
  exp1[is.nan(exp1)] <- 0
  exp2[is.nan(exp2)] <- 0
  
  out <- list(log2((psb_cellprops1+0.0001)/(exp1+0.0001)), log2((psb_cellprops2+0.0001)/(exp2+0.0001)))
  names(out) <- .levels
  
  return(out)
}

observed <- generate_cellcell_profiles(.adj = A, .cellmeta = META, .ptmeta = ptkey, .factor = "line1orf1_group", .levels = c("high", "low"), .randomize = F)
simulated <- replicate(n = 1000, 
                       expr = generate_cellcell_profiles(.adj = A, .cellmeta = META, .ptmeta = ptkey, .factor = "line1orf1_group", .levels = c("high", "low"), .randomize = T), 
                       simplify = F)
simulated <- map(.x = simulated, .f = abind::abind, along = 3)
simulated <- abind::abind(simulated, along = 4)
simulated_diff <- simulated[,,"high",]-simulated[,,"low",]
dimnames(simulated_diff)[[3]] <- 1:1000

get_empircial_p <- function(.index, .neighbor) {
  if (observed_diff[.index, .neighbor] < 0) {
    pval <- mean(simulated_diff[.index, .neighbor, 1:1000] < observed_diff[.index, .neighbor])
  }
  if (observed_diff[.index, .neighbor] > 0) {
    pval <- mean(simulated_diff[.index, .neighbor, 1:1000] > observed_diff[.index, .neighbor])
  }
  return(pval)
}

observed_diff <- (observed$high-observed$low)
pvals <- matrix(data = NA, nrow = 19, ncol = 19, dimnames = list(rownames(observed_diff), colnames(observed_diff)))
for (rn in rownames(pvals)) {
  for (cn in colnames(pvals)) {
    pvals[rn, cn] <- get_empircial_p(rn, cn)
  }
}
pvals

o <- hclust(dist(observed$low))$order
splitter <- rep(c("immune_oi", "not_oi"), times = c(11, 8))
Heatmap(matrix = observed$high[o,o], cluster_rows = F, cluster_columns = F,
        col = colorRamp2(breaks = c(-10, -4, 0, 4, 10), colors = c("blue", scales::muted("blue"), "white", scales::muted("red"), "red")),
        width = ncol(observed$high)*unit(5, "mm"),
        height = nrow(observed$high)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1), 
        row_title = "Index", column_title = "Neighbor", 
        name = "log2ER", 
        column_split = splitter, split = splitter
        ) |> draw(column_title = "LINE1-High")
Heatmap(matrix = observed$low[o,o], cluster_rows = F, cluster_columns = F,
        col = colorRamp2(breaks = c(-10, -4, 0, 4, 10), colors = c("blue", scales::muted("blue"), "white", scales::muted("red"), "red")),
        width = ncol(observed$low)*unit(5, "mm"),
        height = nrow(observed$low)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1), 
        row_title = "Index", column_title = "Neighbor", 
        name = "log2ER", 
        column_split = splitter, split = splitter
        ) |> draw(column_title = "LINE1-Low")

Heatmap(matrix = observed_diff[o,o], cluster_rows = F, cluster_columns = F,
        col = colorRamp2(breaks = c(-10, -4, 0, 4, 10), colors = c("blue", scales::muted("blue"), "white", scales::muted("red"), "red")),
        width = ncol(observed_diff)*unit(5, "mm"),
        height = nrow(observed_diff)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1), 
        row_title = "Index", column_title = "Neighbor", 
        name = "log2FC", 
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(pvals[o,o][i, j] < 0.05) {
            if(observed_diff[o,o][i,j] < -3) {
              grid.text("*", x, y, gp = gpar(col = "white"))
            }
            else {
              grid.text("*", x, y, gp = gpar(col = "black"))
            }
          } 
          else {
            grid.text("", x, y)
          }
        }, 
        column_split = splitter, split = splitter
        )|> draw(column_title = "LINE1-High minus LINE1-Low")

# Now for the overall co-localization score of immune cells
# I am simply getting a sum of the boxes in the larger immune-immune box
immune_coloco_score <- sum(observed_diff[o,o][splitter == "immune_oi", splitter == "immune_oi"])

# I calculate this for each simulation
simulated_scores <- simulated_diff[o,o,][splitter == "immune_oi", splitter == "immune_oi",] |> apply(MARGIN = 3, FUN = sum)

# A really small value means more immune-immune interaction in the low group
ggplot(tibble("score" = simulated_scores)) + 
  geom_histogram(mapping = aes(x = score), bins = 25) + 
  geom_vline(xintercept = immune_coloco_score, colour = "red")

# Here is our p-value: p < 0.1
mean(simulated_scores < immune_coloco_score)

ggplot(tibble("score" = simulated_scores)) + 
  geom_histogram(mapping = aes(x = score), bins = 25) + 
  geom_vline(xintercept = immune_coloco_score, colour = "red") + 
  annotate(geom = "text", x = 100, y = 75, label = paste("p=", mean(simulated_scores < immune_coloco_score), sep = ""))

# Saving what we need for figures
ne_results <- list(
  "obs"=observed,
  "sim"=simulated,
  "pval"=pvals,
  "notes"="index = rows, neighbor = columns"
)
saveRDS(object = ne_results, file = "6-spatial-analysis/line1_patient_level_ne_analysis_results_final_testing.RDS")

# Martin thought that combining some of the immune cell types may be useful as well
cellgroup_map <- c(
  "B-cell"="B/Plasmablast",
  "cholangiocyte"="cholangiocyte", 
  "endothelial"="stromal", 
  "erythrocyte"="erythrocyte", 
  "fibroblast"="stromal",
  "macrophage"="myeloid",
  "mast"="myeloid",
  "mDC"="dendritic",
  "monocyte"="myeloid",
  "neutrophil"="myeloid",
  "NK"="T/NK",
  "pDC"="dendritic",
  "plasmablast"="B/Plasmablast",
  "T CD4 memory"="T/NK",
  "T CD4 naive"="T/NK",
  "T CD8 memory"="T/NK",
  "T CD8 naive"="T/NK",
  "Treg"="T/NK",
  "tumor"="tumor"
)

META$celltype <- plyr::mapvalues(x = META$celltype, from = names(cellgroup_map), to = cellgroup_map)
META$celltype %<>% as.character() %<>% as.factor()

# Doing the analysis again
observed <- generate_cellcell_profiles(.adj = A, .cellmeta = META, .ptmeta = ptkey, .factor = "line1orf1_group", .levels = c("high", "low"), .randomize = F)
simulated <- replicate(n = 1000, 
                       expr = generate_cellcell_profiles(.adj = A, .cellmeta = META, .ptmeta = ptkey, .factor = "line1orf1_group", .levels = c("high", "low"), .randomize = T), 
                       simplify = F)
simulated <- map(.x = simulated, .f = abind::abind, along = 3)
simulated <- abind::abind(simulated, along = 4)
simulated_diff <- simulated[,,"high",]-simulated[,,"low",]
dimnames(simulated_diff)[[3]] <- 1:1000

get_empircial_p <- function(.index, .neighbor) {
  if (observed_diff[.index, .neighbor] < 0) {
    pval <- mean(simulated_diff[.index, .neighbor, 1:1000] < observed_diff[.index, .neighbor])
  }
  if (observed_diff[.index, .neighbor] > 0) {
    pval <- mean(simulated_diff[.index, .neighbor, 1:1000] > observed_diff[.index, .neighbor])
  }
  return(pval)
}

observed_diff <- (observed$high-observed$low)
pvals <- matrix(data = NA, nrow = 8, ncol = 8, dimnames = list(rownames(observed_diff), colnames(observed_diff)))
for (rn in rownames(pvals)) {
  for (cn in colnames(pvals)) {
    pvals[rn, cn] <- get_empircial_p(rn, cn)
  }
}
pvals

o <- hclust(dist(observed$low))$order
Heatmap(matrix = observed$high[o,o], cluster_rows = F, cluster_columns = F,
        col = colorRamp2(breaks = c(-10, -4, 0, 4, 10), colors = c("blue", scales::muted("blue"), "white", scales::muted("red"), "red")),
        width = ncol(observed$high)*unit(5, "mm"),
        height = nrow(observed$high)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1), 
        row_title = "Index", column_title = "Neighbor", 
        name = "log2ER"
) |> draw(column_title = "LINE1-High")
Heatmap(matrix = observed$low[o,o], cluster_rows = F, cluster_columns = F,
        col = colorRamp2(breaks = c(-10, -4, 0, 4, 10), colors = c("blue", scales::muted("blue"), "white", scales::muted("red"), "red")),
        width = ncol(observed$low)*unit(5, "mm"),
        height = nrow(observed$low)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1), 
        row_title = "Index", column_title = "Neighbor", 
        name = "log2ER", 
) |> draw(column_title = "LINE1-Low")

Heatmap(matrix = observed_diff[o,o], cluster_rows = F, cluster_columns = F,
        col = colorRamp2(breaks = c(-10, -4, 0, 4, 10), colors = c("blue", scales::muted("blue"), "white", scales::muted("red"), "red")),
        width = ncol(observed_diff)*unit(5, "mm"),
        height = nrow(observed_diff)*unit(5, "mm"), 
        rect_gp = gpar(col = "grey", lwd = 1), 
        row_title = "Index", column_title = "Neighbor", 
        name = "log2FC", 
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(pvals[o,o][i, j] < 0.05) {
            if(observed_diff[o,o][i,j] < -3) {
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
)|> draw(column_title = "LINE1-High minus LINE1-Low")

# This is interesting: B/Plasma -- B/Plasma co-localization is enriched in the low group

# Saving what we need for figures
ne_results_collapsed <- list(
  "obs"=observed,
  "sim"=simulated,
  "pval"=pvals,
  "notes"="index = rows, neighbor = columns"
)
saveRDS(object = ne_results, file = "6-spatial-analysis/line1_patient_level_ne_analysis_results_final_testing_collapsed.RDS")

ssseeeddd
# [1]       10403         624   168352879  1488529380 -1691384427   879462738  -110616405   582342640  -544924879

