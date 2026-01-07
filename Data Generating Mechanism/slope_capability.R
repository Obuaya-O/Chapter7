# Slope capacity

# Aim: Finding the biomarkers and outcomes with two distinct time points for a slope calculation

# Library
library(tidyverse)

#Dataset
master_data <- read_csv("adni_master_longitudinal_cohort.csv")

#Biomarkers
biomarkers <- c("MMSCORE", "TOTAL13", "ABETA40", "ABETA42","PTAU", "CENTILOIDS",
                "COMPOSITE_REF_VOLUME", "TAU_META_ROI", "HIPPOCAMPUS_VOL",
                "AB42_F", "AB40_F", "PLASMAPTAU181")

# Count how many subjects have >= 2 points
check_overlap <- function(df, var1, var2) {
  # Error handling if columns don't exist
  if (!var1 %in% names(df)) return(paste(var1, "not found"))
  if (!var2 %in% names(df)) return(paste(var2, "not found"))
  
  count <- df %>%
    group_by(RID) %>%
    summarise(
      has_v1 = sum(!is.na(.data[[var1]])) >= 2,
      has_v2 = sum(!is.na(.data[[var2]])) >= 2
    ) %>%
    filter(has_v1 & has_v2) %>%
    nrow()
  
  return(count)
}

# Check overlap for D-matrix
cat("\n--- OVERLAP CHECK: ADAS-Cog vs. Biomarkers ---\n")
for (bio in biomarkers) {
  n <- check_overlap(master_data, "TOTAL13", bio)
  cat("TOTAL13 <->", bio, ":\t", n, "subjects\n")
}

cat("\n--- OVERLAP CHECK: MMSE vs. Biomarkers ---\n")
for (bio in biomarkers) {
  n <- check_overlap(master_data, "MMSCORE", bio)
  cat("MMSCORE <->", bio, ":\t", n, "subjects\n")
}

cat("\n--- OVERLAP CHECK: Cognitive Anchor ---\n")
n_cog <- check_overlap(master_data, "TOTAL13", "MMSCORE")
cat("TOTAL13 <-> MMSCORE :\t", n_cog, "subjects\n")
