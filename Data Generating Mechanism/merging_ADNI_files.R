# Merging ADNI Datasets to Create Master Longitudinal Cohort
# Data is available at: https://adni.loni.usc.edu/

# Packages
library(tidyverse)

# Load main dataset
adni_inter <- read_csv(".csv") # Interim dataset without biomarker data. Due to data sharing restrictions, this file is not provided here.

# Load biomarker datasets
ab_pet <- read_csv(".csv")
tau_pet <- read_csv(".csv")
ab_plasma <- read_csv(".csv")
ptau_plasma <- read_csv(".csv")

# Function to clean keys AND remove duplicates
clean_and_deduplicate <- function(df) {
  df %>%
    mutate(
      RID = as.numeric(RID),
      VISCODE = str_trim(as.character(VISCODE))
    ) %>%
    distinct(RID, VISCODE, .keep_all = TRUE)
}

# Clean the interim cohort
adni_inter_clean <- adni_inter %>%
  mutate(
    RID = as.numeric(str_sub(subject_id, -4)) 
  ) %>%
  rename(VISCODE = visit) %>% #Keep consistent naming
  clean_and_deduplicate() # Remove duplicates 

# Clean and de-duplicate all biomarker files
ab_pet_clean    <- clean_and_deduplicate(ab_pet)
tau_pet_clean   <- clean_and_deduplicate(tau_pet)
ab_plasma_clean <- clean_and_deduplicate(ab_plasma)
ptau_plasma_clean <- clean_and_deduplicate(ptau_plasma)

# Merging the datasets
master_longitudinal_cohort <- adni_inter_clean %>% 
  left_join(select(ab_pet_clean, RID, VISCODE, CENTILOIDS, AMYLOID_STATUS,
                   COMPOSITE_REF_VOLUME),
            by = c("RID", "VISCODE")) %>% 
  left_join(select(tau_pet_clean, RID, VISCODE, TAU_META_ROI = META_TEMPORAL_SUVR, HIPPOCAMPUS_VOL = HIPPOCAMPUS_VOLUME),
            by = c("RID", "VISCODE")) %>% 
  left_join(select(ab_plasma_clean, RID, VISCODE, AB42_F, AB40_F),
            by = c("RID", "VISCODE")) %>% 
  left_join(select(ptau_plasma_clean, RID, VISCODE, PLASMAPTAU181),
            by = c("RID", "VISCODE"))

# Check the results
cat("Original Subjects:", n_distinct(adni_inter_clean$RID), "\n")
cat("Final Subjects:   ", n_distinct(master_longitudinal_cohort$RID), "\n")
cat("Final Rows:       ", nrow(master_longitudinal_cohort), "\n")

# Save dataset
write_csv(master_longitudinal_cohort, "adni_master_longitudinal_cohort.csv")

