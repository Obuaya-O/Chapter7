# setup_params.R
# PURPOSE: Define the 81 Scenarios for the Thesis Simulation
# 1. Sample Size (3 levels)
# 2. Effect Size (Null, Medium=25%, Large=50%)
# 3. Correlation Level (Real, Moderate=0.4, Strong=0.6)
# 4. Missing Data (MCAR, MAR, MNAR)

df <- expand.grid(
  sample_size  = c(260, 530, 1072),
  effect_size  = c("Null", "Medium", "Large"),
  corr_level   = c("Real", "Moderate", "Strong"),
  missing_type = c("MCAR", "MAR", "MNAR"),
  chunk_id     = 1:20 
)

df$scenario_id <- 1:nrow(df)
saveRDS(df, "simulation_parameters.rds")

print(paste("Created simulation_parameters.rds with", nrow(df), "scenarios."))
print(head(df))
