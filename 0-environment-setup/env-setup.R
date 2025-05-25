# HCC TMA Project Environment Setup
# This should work on an ARM mac

renv::init()
.libPaths()

# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project-final/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"

# Should work for any device
renv::install(packages = "BiocManager")
renv::install(packages = "bioc::SingleCellExperiment")
renv::install(packages = "lme4")
renv::install(packages = "lmerTest")

# May need to install a different binary that is suited for your device
install.packages("https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.4/Rfast_2.1.3.tgz", repos = NULL, type = "binary")

# Should work for any device
renv::install(packages = "NanoString-BioStats/InSituCor")
renv::install(packages = "Seurat")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
renv::install("nebula")
renv::install("Bioc::ComplexHeatmap")
renv::install("immunogenomics/presto@glmm")
renv::install("immunogenomics/singlecellmethods")
renv::install("Bioc::limma")
renv::install("Bioc::edgeR")
renv::install("Bioc::DESeq2")
renv::install("Bioc::glmGamPoi")
renv::install(packages = list("ggdendro", "ggh4x"))
renv::install("bioc::sparseMatrixStats")
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
install.packages("remotes")

# May depend on your compiler. I used the homebrew gcc, evidenced here: 
system("which gfortran")
# /opt/homebrew/bin/gfortran

# Therefore, I did this: 
dir.create('~/.R')
file.create('~/.R/Makevars')
# Opened the Makevars file and added these lines to it (I have version 14.2.0_1): 
# FC = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
# F77 = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
# FLIBS = -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14
renv::install("NanoString-Biostats/InSituType")

# See https://mac.r-project.org/tools/ for gcc install not from homebrew. 
# I believe that you can install it and specify the corresponding paths in the Makevars file.

# Should work, if you have the package in the same directory as this script.
install.packages("0-environment-setup/smiDE_0.0.2.01.tar.gz", repos = NULL, type = "source")

# Should work on any device
renv::install("ggpubr")
renv::install("jpeg")
renv::install("openxlsx")
renv::install("arrow")
renv::install("davidsjoberg/ggsankey")
renv::install("bioc::biomaRt")
renv::install("R.utils")
renv::install("emmeans")

# If all of these packages are installed, you should not experience any problems.
renv::snapshot()

# Restart R as a precaution.
.rs.restartR(clean = T)

