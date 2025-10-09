# Library
library(dplyr)

# Load datasets
main_adni <- read_csv("adni_matched_cohort 2_current.csv")
image_adni<- read_csv("adni_imaging_ds.csv")

# Columns to keep
keep_main <- c(
  "subject_id", "visit", "PTGENDER", "GENOTYPE", "entry_age", "PTEDUCAT", "MMSCORE",
  "KEYMED", "TOTAL13", "ABETA40", "ABETA42", "PTAU", "PLASMA_BIOMARKER", "TESTVALUE"
)

keep_image <- c(
  "PTID", "VISCODE", "COMPOSITE_REF_VOLUME" 
)

#Apply the selection
main_adni_clean <- main_adni %>% select(all_of(keep_main))
image_adni_clean <- image_adni %>% select(all_of(keep_image))

#Merge the tables
image_adni_clean <- image_adni_clean %>% 
  rename(
    subject_id = PTID,
    visit = VISCODE
  )

adni_matched_dataset <- left_join(main_adni_clean, image_adni_clean, by = c("subject_id", "visit"))
# Have a look-see
glimpse(adni_matched_dataset)

#save it
write.csv(adni_matched_dataset, "adni_matched_dataset.csv")
