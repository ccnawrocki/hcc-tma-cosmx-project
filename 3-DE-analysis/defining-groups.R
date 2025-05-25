# Group Definitions #
## Cole Nawrocki ##

# Environment
rm(list = ls())
.libPaths()

# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project-final/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"    

# Packages
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggprism)
library(circlize)
library(Seurat)

# Data
tmas <- LoadSeuratRds(file = "hcc_tmas_final.RDS")
clindb <- openxlsx::read.xlsx("~/Partners HealthCare Dropbox/Ting lab/Mouse metastasis Repeats Avril Coley/HCC Repeat Paper/HCC Repeat Database 2025.xlsx", sheet = 1, rows = 1:40, cols = 2:47)

# Making core to path id key
keys <- list()
for (tma_nm in c("tma1", "tma2")) {
  keys[[tma_nm]] <- openxlsx::read.xlsx("1-QC-and-processing/hcc_tma_maps.xlsx", sheet = paste(tma_nm, "map", sep = "-"), cols = 11:12)
}
pathkey <- bind_rows(keys, .id = "slide")
pathkey$core_global <- paste(stringr::str_sub(string = pathkey$slide, start = 4, end = 4), pathkey$core, sep = "")
pathkey$patient <- stringr::str_split(string = pathkey$pathology, pattern = " ", simplify = T)[,1]

# Subset down to only the cores that are in the cosmx data
tmas$core_global <- paste(stringr::str_sub(string = tmas$slide, start = 8, end = 8), tmas$core, sep = "")
incosmx <- (tmas$core_global |> unique())

# Add MRNs and de-identified patient IDs
pathkey <- filter(pathkey, core_global %in% incosmx)
pathkey$MRN <- plyr::mapvalues(x = pathkey$patient, from = clindb$patient, to = clindb$MRN)
missed <- readRDS(file = "3-DE-analysis/missed_mrns.RDS") # Saved these here, so this script has no PII in it
pathkey[pathkey$patient %in% names(missed),]$MRN <- missed[pathkey[pathkey$patient %in% names(missed),]$patient]
pathkey$patient_deid <- plyr::mapvalues(x = pathkey$core_global, from = tmas$core_global, to = as.character(tmas$patient))

# Double checking that our de-identified patient IDs are okay... they are.
pathkey |> group_by(MRN) |> tally() |> pull(n)
# [1] 2 1 1 1 1 1 1 2 1 3 1 1 1 1 1 1 1 2 3 1 1 1 1
pathkey |> group_by(MRN, patient_deid) |> tally() |> pull(n)
# [1] 2 1 1 1 1 1 1 2 1 3 1 1 1 1 1 1 1 2 3 1 1 1 1

# LINE1 tumor counts per million
d1 <- tmas@meta.data |>
  group_by(core_global, celltype) |>
  filter(celltype == "tumor") |>
  dplyr::summarise(patient = unique(patient) |> as.character(),
                   slide = unique(slide) |> as.character(),
                   total_tumor_counts = sum(nCount_RNA))
pathkey$total_tumor_counts <- plyr::mapvalues(x = pathkey$core_global, from = d1$core_global, to = d1$total_tumor_counts) |> as.numeric()

d2 <- tmas@assays$RNA@counts["LINE1-ORF1", (tmas$celltype == "tumor"), drop = F]
mm <- model.matrix(~0+core_global, data = tmas@meta.data[(tmas$celltype == "tumor"),])
l1cts <- (d2 %*% mm) |> as.matrix() |> t() |> as.data.frame()
l1cts$core_global <- gsub(x = rownames(l1cts), pattern = "core_global", replacement = "")
pathkey$total_tumor_line1_counts <- plyr::mapvalues(x = pathkey$core_global, from = l1cts$core_global, to = l1cts$`LINE1-ORF1`) |> as.numeric()

ptkey <- group_by(pathkey, MRN) |> 
  summarise(patient = unique(patient), 
            patient_deid = unique(patient_deid), 
            total_tumor_counts = sum(total_tumor_counts), 
            total_tumor_line1_counts = sum(total_tumor_line1_counts))
ptkey$line1orf1_tumor_cpm <- (ptkey$total_tumor_line1_counts/ptkey$total_tumor_counts)*1e6 

# Plotting
ptkey$patient_deid <- factor(x = ptkey$patient_deid, levels = ptkey |> arrange(desc(line1orf1_tumor_cpm)) |> pull(patient_deid))
ggplot() + 
  geom_bar(data = ptkey, mapping = aes(x = patient_deid, y = line1orf1_tumor_cpm), stat = "identity") + 
  geom_hline(yintercept = quantile(x = ptkey$line1orf1_tumor_cpm, probs = 2/3), linetype = "dashed", color = "red")

# Defining groups by core
ptkey$line1orf1_group <- ifelse(test = ptkey$line1orf1_tumor_cpm > quantile(x = ptkey$line1orf1_tumor_cpm, probs = 2/3), 
                                yes = "high", no = "low")
ggplot() + 
  geom_bar(data = ptkey, mapping = aes(x = patient_deid, y = line1orf1_tumor_cpm, fill = line1orf1_group), stat = "identity") + 
  geom_hline(yintercept = quantile(x = ptkey$line1orf1_tumor_cpm, probs = 2/3), linetype = "dashed", color = "red") + 
  scale_fill_manual(values = c("high" = "dodgerblue", "low" = "lightblue")) + 
  theme_bw()

# Notes: 
## There are 6 patients in the CosMx data who exist in the clinical database but who have no LINE1 staining, likely due to lack of tissue.
## One of these patients passed away 15 days after surgery.
## So, the clinical cohort and the CosMx cohort overlap, but not very much.

# Saving what we may need downstream
## CosMx patients
fulldb <- openxlsx::read.xlsx("~/Partners HealthCare Dropbox/Cole Nawrocki/HCC Repeat Paper/HCC Repeat Database 2025.xlsx", sheet = "Database")

# Differentiation
ptkey$differentiation <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$`Differentiation.Well,.Moderate-Well,.moderate,.Moderate-Poor,.Poor)`)
ptkey[grepl(pattern = "^[0-9]", x =  ptkey$differentiation),]$differentiation <- 
  plyr::mapvalues(x = ptkey$differentiation[grepl(pattern = "^[0-9]", x =  ptkey$differentiation)], from = fulldb$MRN, to = fulldb$differentiation)

# Five year OS
ptkey$five_year_post_surgery_OS <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$os_days_fromsurgery_5yearsmax) |> as.numeric()
ptkey[ptkey$five_year_post_surgery_OS > 10000,]$five_year_post_surgery_OS <- 
  plyr::mapvalues(x = ptkey$five_year_post_surgery_OS[ptkey$five_year_post_surgery_OS > 10000], from = fulldb$MRN, to = fulldb$os_days_fromsurgery_5yearsmax) |> as.numeric()

# Death and censor
ptkey$death <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$death1censored05years) |> as.numeric()
ptkey[ptkey$death > 10,]$death <- 
  plyr::mapvalues(x = ptkey$death[ptkey$death > 10], from = fulldb$MRN, to = fulldb$death1.censored0.5years) |> as.numeric()

# Primary etiology
ptkey$primary_etiology <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$prim_etiology1)
ptkey[grepl(pattern = "^[0-9]", x =  ptkey$primary_etiology),]$primary_etiology <- 
  plyr::mapvalues(x = ptkey$primary_etiology[grepl(pattern = "^[0-9]", x =  ptkey$primary_etiology)], from = fulldb$MRN, to = fulldb$etiology1)
ptkey$primary_etiology <- ifelse(test = ptkey$primary_etiology == ".", yes = NA, no = ptkey$primary_etiology)

# Surgery
ptkey$surgery <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$`Surgery.(Resection,.Transplant)`)
ptkey[grepl(pattern = "^[0-9]", x =  ptkey$surgery),]$surgery <- 
  plyr::mapvalues(x = ptkey$surgery[grepl(pattern = "^[0-9]", x =  ptkey$surgery)], from = fulldb$MRN, to = fulldb$sample)

# Sex
ptkey$sex <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$Sex_1M_2F) |> as.numeric()
ptkey[ptkey$sex > 10,]$sex <- 
  plyr::mapvalues(x = ptkey$sex[ptkey$sex > 10], from = fulldb$MRN, to = fulldb$Gender_1M_2F) |> as.numeric()
ptkey$sex <- ifelse(test = ptkey$sex == 1, yes = "male", no = "female")

# Age
ptkey$age_at_surgery <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$Age) |> as.numeric()
ptkey[ptkey$age_at_surgery > 10000,]$age_at_surgery <- 
  plyr::mapvalues(x = ptkey$age_at_surgery[ptkey$age_at_surgery > 10000], from = fulldb$MRN, to = fulldb$Age_surgery) |> as.numeric()

# Tumor size 
ptkey$tumor_size_cm <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$`Tumor.size.(cm,.median,.range)`) |> as.numeric()
ptkey[ptkey$tumor_size_cm > 10000,]$tumor_size_cm <- 
  plyr::mapvalues(x = ptkey$tumor_size_cm[ptkey$tumor_size_cm > 10000], from = fulldb$MRN, to = fulldb$tumor_size) |> as.numeric()

# Serum AFP
ptkey$serum_afp_ng_per_mL <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$`Serum.AFP.ng/mL.(Median,.Range)`) |> as.numeric()
ptkey[ptkey$serum_afp_ng_per_mL > 100000,]$serum_afp_ng_per_mL <- 
  plyr::mapvalues(x = ptkey$serum_afp_ng_per_mL[ptkey$serum_afp_ng_per_mL > 100000], from = fulldb$MRN, to = fulldb$serum_afp) |> as.numeric()

# tstage
ptkey$t_stage <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$tstage) |> as.numeric()
ptkey[ptkey$t_stage > 10000,]$t_stage <- 
  plyr::mapvalues(x = ptkey$t_stage[ptkey$t_stage > 10000], from = fulldb$MRN, to = fulldb$tstage) |> as.numeric()

# in ISH data?
ptkey$in_ish_data <- (ptkey$MRN %in% clindb$MRN)
ptkey$ish_line1_counts_per_square_um <- plyr::mapvalues(x = ptkey$MRN, from = clindb$MRN, to = clindb$LINE1_count_density_sq_um) |> as.numeric()
ptkey$ish_line1_counts_per_square_um[ptkey$MRN == ptkey$ish_line1_counts_per_square_um] <- NA

# Saving
write.csv(x = ptkey, file = "3-DE-analysis/cosmx-patient-data.csv")
write.csv(x = ptkey |> select(-MRN, -patient), file = "3-DE-analysis/cosmx-patient-data-de-identified.csv")


## ISH patients
ishpts <- clindb |> 
  select(1, 2, 5, 6, 7, 8, 9, 10, 15, 16, 19, 21, 22, 23, 24, 25, 31:39, 41, 42, 43) 
colnames(ishpts) <- c("patient", "MRN", "tma", "core", "line1_counts_per_square_um", "line1_area_per_square_um", 
                      "five_year_post_surgery_OS", "death", 
                      "sex", "age_surgery", "primary_etiology", "surgery", "tumor_size", "differentiation", "serum_afp_ng_per_mL", "tstage",
                      "CD3_10000", "CD4_10000", "CD8_10000", "PD1_10000", "LAG3_10000", "TIM3_10000", "FOXP3_10000", "CD163_10000", "CD20_10000", 
                      "HERVK_PerCell", "HSATII_PerCell", "HERVH_PerCell")
ishpts$sex <- ifelse(test = ishpts$sex == 1, yes = "male", no = "female")
ishpts$in_cosmx_data <- (ishpts$MRN %in% ptkey$MRN)
ishpts$cosmx_line1orf1_tumor_cpm <- plyr::mapvalues(x = ishpts$MRN, from = ptkey$MRN, to = ptkey$line1orf1_tumor_cpm) |> as.numeric()
ishpts$cosmx_line1orf1_tumor_cpm[ishpts$MRN == ishpts$cosmx_line1orf1_tumor_cpm] <- NA

# Saving 
ishpts$patient_deid <- plyr::mapvalues(x = ishpts$MRN, from = ptkey$MRN, to = as.character(ptkey$patient_deid))
ishpts <- ishpts |> split(f = ishpts$in_cosmx_data)
ishpts$`FALSE`$patient_deid %<>% as.factor() %<>% as.numeric()
ishpts$`FALSE`$patient_deid <- paste(ishpts$`FALSE`$tma, ishpts$`FALSE`$patient_deid, sep = "p")
ishpts <- bind_rows(ishpts)
write.csv(x = ishpts, file = "3-DE-analysis/ish-patient-data.csv")
write.csv(x = ishpts |> select(-MRN, -patient), file = "3-DE-analysis/ish-patient-data-de-identified.csv")

## Info for Figures ------------------------------------------------------------
# CosMx N patients: 23
ptkey |> nrow()

# Overlap with ISH data N: 17
intersect(ptkey$MRN, ishpts$MRN) |> n_distinct()

# Patients added to ISH data N: 22
nrow(ishpts)-17

# So the paper should convey that 48 cores were embedded in FFPE on TMA. After setting up the slides for CosMx, 
# removing controls, and QC filtering, 30 cores remained, representing 23 patients. After the CosMx analysis, 
# ISH and IHC data was collected after serial sectioning for as many patients as possible. 17 of the patients 
# included in the CosMx cohort were able to be included in ISH staining. An additional 22 patients were included 
# in the ISH staining cohort. So, overall: N=23 for CosMx and N=39 for ISH. 17 patients are included in both 
# datasets. We will need summary tables for both datasets.



