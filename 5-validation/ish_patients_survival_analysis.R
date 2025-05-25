# Packages
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggsurvfit)
library(survival)
library(survminer)

pdf(file = "5-clinical-analysis/ish_patients_survival_analysis.pdf", width = 6, height = 6)

# Data
ishpts <- read.csv("3-DE-analysis/ish-patient-data.csv", row.names = 1)

# Groups: Upper tercile vs. bottom two terciles
ishpts$patient <- factor(x = ishpts$patient, levels = ishpts |> arrange(desc(line1_counts_per_square_um)) |> pull(patient))
ishpts$line1_group <- ifelse(test = ishpts$line1_counts_per_square_um > quantile(x = ishpts$line1_counts_per_square_um, probs = 2/3), 
                             yes = "high", no = "low")
ggplot() + 
  geom_bar(data = ishpts, mapping = aes(x = patient, y = line1_counts_per_square_um, fill = line1_group), stat = "identity") + 
  geom_hline(yintercept = quantile(x = ishpts$line1_counts_per_square_um, probs = 2/3), linetype = "dashed", color = "red") + 
  scale_fill_manual(values = c("high" = "firebrick", "low" = "dodgerblue")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 0))

# Fitting the survival model
fit <- survfit(Surv(five_year_post_surgery_OS, death) ~ line1_group,
               data = ishpts)

# Simple plot
ggsurvplot(fit = fit, data = ishpts, risk.table = T, conf.int = F, pval = T, palette = c("firebrick", "dodgerblue"))

# More customizable plot
fit2 <- survfit2(Surv(five_year_post_surgery_OS, death) ~ line1_group,
                 data = ishpts) 
ggsurvfit(x = fit2, linewidth = 1) +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) +
  #add_risktable() + 
  add_pvalue(location  = "annotation") +
  scale_color_manual(values = c("high" = "firebrick", "low" = "dodgerblue")) + 
  theme_prism()

dev.off()

