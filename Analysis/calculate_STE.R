# Library
library(tidyverse)

# 1. Load your combined data 

# 2. Function to calculate STE from a linear model
calculate_ste_root <- function(df, outcome_type) {
  
  # Fit Model: Clinical Effect ~ Surrogate Effect
  # We use the 'multi' surrogate (panel) or 'bench' (amyloid) depending on what data is passed
  model <- lm(est_bench_clin ~ est_multi_surr, data = df)
  
  # We need to find the surrogate value (x) where the Prediction Interval crosses 0.
  # Quadratic equation (prediction interval)-> we solve it numerically.
  
  # Define the prediction function for the Lower/Upper bound
  # ADAS: Improvement is NEGATIVE. We want the Upper Bound of the prediction to be < 0.
  # MMSE: Improvement is POSITIVE. We want the Lower Bound of the prediction to be > 0.
  
  pred_fun <- function(x) {
    p <- predict(model, newdata = data.frame(est_multi_surr = x), interval = "prediction", level = 0.95)
    if (outcome_type == "ADAS") return(p[1, "upr"]) # Upper bound must be < 0
    if (outcome_type == "MMSE") return(p[1, "lwr"]) # Lower bound must be > 0
  }
  
  # Search for the root (where bound = 0)
  # We search a reasonable range of surrogate effects (e.g., -5 to 5 SDs or Centiloids)
  tryCatch({
    # Attempt to find root
    if (outcome_type == "ADAS") {
      # For ADAS, we expect negative surrogate values to predict clinical benefit
      root <- uniroot(pred_fun, c(-10, 0))$root 
    } else {
      # For MMSE, we expect positive surrogate values (or negative if using Amyloid benchmarks)
      # Assuming Composite is aligned with MMSE (Higher = Better), range is 0 to 10
      root <- uniroot(pred_fun, c(0, 10))$root
    }
    return(root)
  }, error = function(e) {
    return(NA) # Returns NA if STE is infinite or outside range (common in poor surrogates)
  })
}

# 3. Apply Calculation
# We group by Model, N, Corr, Miss.
# IMPORTANT: We use ALL effect sizes (Null, Medium, Large) together to fit the regression.
ste_results <- all_data %>%
  group_by(Model, N, Corr, Miss) %>%
  summarise(
    # Check if Outcome is ADAS or MMSE for directionality
    Outcome_Type = ifelse(grepl("ADAS", unique(Model)), "ADAS", "MMSE"),
    
    # Calculate STE for the multibiomarker panel
    STE_Panel = calculate_ste_root(cur_data(), Outcome_Type),
    
    # Calculate STE for the benchmark (amyloid)
    # Note: We need to temporarily rename columns to use the same function
    STE_Bench = {
      df_bench <- cur_data()
      df_bench$est_multi_surr <- df_bench$est_bench_surr # Swap variable for function
      calculate_ste_root(df_bench, Outcome_Type)
    },
    .groups = "drop"
  )

# 4. View and Save
print(head(ste_results))
write.csv(ste_results, "Table_STE_Final.csv", row.names = FALSE)

# Interpretation:
# The STE value is the unit change in the surrogate needed to guarantee clinical success.
# For weighted panel (standardised), an STE of -1.5 means "You need to shift the score by -1.5 SDs".
# For benchmark (amyloid), an STE of -20 means "You need to lower amyloid by 20 Centiloids".