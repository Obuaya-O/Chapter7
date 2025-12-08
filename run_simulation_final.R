#!/usr/bin/env Rscript

# PURPOSE: Run the "In Silico" Phase III Trial Simulation
# INPUT:   simulation_parameters.rds
# OUTPUT:  simulation_results_N1000_jobID.rds

# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(joineRML)
  library(parallel)
  library(MASS) # For mvrnorm
  library(Matrix)
})

# 1. READ COMMAND LINE ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)
# If no argument provided (default to 1 for testing)
job_id <- if (length(args) == 0) 1 else as.integer(args[1])

# 2. LOAD PARAMETERS (CRITICAL FIX)
if (!file.exists("simulation_parameters.rds")) {
  stop("Error: simulation_parameters.rds not found.")
}
params <- readRDS("simulation_parameters.rds")

# 3. SETUP SIMULATION CONSTANTS
n_reps <- 1000  
cat("----------------------------------------------------------\n")
cat("STARTING SIMULATION | Job ID:", job_id, "\n")
cat("Target Replicates:", n_reps, "\n")
cat("----------------------------------------------------------\n")

# 4. DEFINE SCENARIOS
scenarios <- expand.grid(
  sample_size   = c(260, 530, 1072),
  effect_size   = c("Null", "Medium", "Large"),
  missing_type  = c("MCAR", "MAR", "MNAR"),
  corr_strength = c("Weak", "Moderate", "Strong")
) %>% 
  mutate(scenario_id = row_number()) # Add ID here

# Select ONLY the scenario for this specific job
scen <- scenarios[job_id, ] 

print(paste("Running Scenario:", scen$scenario_id))
print(scen)

# 5. RUN PARALLEL SIMULATION
n_cores <- detectCores() - 1

# mclapply returns a LIST of results automatically
results_list <- mclapply(1:n_reps, function(i) {
    
    # --- A. MODIFY D-MATRIX (Sensitivity Analysis) ---
    D_current <- params$D_matrix 
    
    # Define Target Correlations
    target_ab_corr <- case_when(
      scen$corr_strength == "Weak"     ~ -0.24, 
      scen$corr_strength == "Moderate" ~ -0.40,
      scen$corr_strength == "Strong"   ~ -0.60
    )
    
    target_multi_corr <- case_when(
      scen$corr_strength == "Weak"     ~ -0.20,
      scen$corr_strength == "Moderate" ~ -0.45,
      scen$corr_strength == "Strong"   ~ -0.65
    )
    
    # Helper Function
    update_covariance <- function(D, var1, var2, target_r) {
      idx1 <- which(rownames(D) == var1)
      idx2 <- which(rownames(D) == var2)
      if(length(idx1) == 0 || length(idx2) == 0) return(D)
      
      sd1 <- sqrt(D[idx1, idx1])
      sd2 <- sqrt(D[idx2, idx2])
      
      original_cov <- D[idx1, idx2]
      # If original covariance is 0, assume negative correlation (Clinical vs Bio)
      original_sign <- if(sign(original_cov) == 0) -1 else sign(original_cov)
      
      new_cov <- abs(target_r) * original_sign * sd1 * sd2
      
      D[idx1, idx2] <- new_cov
      D[idx2, idx1] <- new_cov
      return(D)
    }
    
    # Apply tweaks
    D_current <- update_covariance(D_current, "ADAS_Cog", "AB_PET", target_ab_corr)
    
    # Ensure these names match your params$D_matrix exactly
    biomarkers <- c("tau_PET", "vMRI", "CSF_ptau", "Plasma_ptau", "CSF_AB42")
    for(bm in biomarkers) {
      D_current <- update_covariance(D_current, "ADAS_Cog", bm, target_multi_corr)
    }
    
    # --- B. GENERATE DATA ---
    n_subj <- scen$sample_size
    
    # Matrix safety check
    eigen_values <- eigen(D_current, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigen_values <= 1e-06)) {
      D_current <- as.matrix(Matrix::nearPD(D_current, corr = FALSE, keepDiag = TRUE)$mat)
    }

    random_effects <- suppressWarnings(
      mvrnorm(n = n_subj, mu = rep(0, ncol(D_current)), Sigma = D_current)
    )
    colnames(random_effects) <- colnames(D_current)
    random_effects <- as.data.frame(random_effects)
    random_effects$sub_id <- 1:n_subj
    
    # Time points
    time_points <- seq(0, 2, by = 0.5)
    
    sim_data <- expand.grid(sub_id = 1:n_subj, time = time_points) %>%
      left_join(random_effects, by = "sub_id")
    
    # Assign Groups
    subj_info <- data.frame(sub_id = 1:n_subj) %>%
      mutate(group = sample(c(0, 1), n_subj, replace = TRUE))
    
    sim_data <- sim_data %>% left_join(subj_info, by = "sub_id")
    
    # Calculate Outcomes
    outcomes <- colnames(D_current)
    
    for(out in outcomes) {
      b0 <- params$intercepts[out]
      b1 <- params$beta[out]
      sigma <- params$sigma[out]
      
      pct_slowing <- case_when(
        scen$effect_size == "Null"   ~ 0.00,
        scen$effect_size == "Medium" ~ 0.25,
        scen$effect_size == "Large"  ~ 0.50
      )
      
      effective_slope <- ifelse(sim_data$group == 1, b1 * (1 - pct_slowing), b1)
      rand_slope_col <- sim_data[[out]]
      noise <- rnorm(nrow(sim_data), mean=0, sd=sigma)
      
      sim_data[[out]] <- b0 + (effective_slope * sim_data$time) + 
        (rand_slope_col * sim_data$time) + noise
    }
    
    # C. MISSINGNESS 
    if(scen$missing_type != "MCAR") { # Assuming MCAR logic for simplicity here
       # Simple random dropout
       sim_data <- sim_data %>%
         mutate(
           is_missing = runif(n()) < 0.05,
           ADAS_Cog = ifelse(is_missing & time > 0, NA, ADAS_Cog)
         ) %>%
         dplyr::select(-is_missing)
    }

    #  D. CREATE COMPOSITE SCORE 
    sim_data <- sim_data %>%
      mutate(
        z_Tau    = as.numeric(scale(tau_PET)),
        z_vMRI   = as.numeric(scale(vMRI)),
        z_Plasma = as.numeric(scale(Plasma_ptau)),
        z_csfAB  = as.numeric(scale(CSF_AB42)),
        z_csfptau= as.numeric(scale(CSF_ptau)),
        
        # Note: vMRI (thickness) goes DOWN as disease progresses. 
        # So low vMRI = Bad. High Z = Good. 
        # We want High Score = Bad. So we invert vMRI (-z_vMRI).
        Score_Multi = (0.56 * z_Tau) + (0.55 * -z_vMRI) + (0.40 * z_Plasma) + 
                      (0.28 * -z_csfAB) + (0.04 * z_csfptau)
      )
    
    # E. FIT MODELS (mjoin) 
    # FIX: Function is mjoin(), not joineRML()
    # ADAS-Cog models

    # Model 1: Benchmark
    fit_bench <- tryCatch({
      mjoin(
        data = sim_data,
        Y.names = c("ADAS_Cog", "AB_PET"), 
        time = "time", 
        long.formulas = list(
          ADAS_Cog ~ 1 + time * group,
          AB_PET   ~ 1 + time * group
        ),
        d = list(~1+time, ~1+time),
        control = list(allow.n.na = TRUE) # Allow missing data
      )
    }, error = function(e) return(NULL))
    
    # Model 2: Multibiomarker
    fit_multi <- tryCatch({
      mjoin(
        data = sim_data,
        Y.names = c("ADAS_Cog", "Score_Multi"), 
        time = "time",
        long.formulas = list(
          ADAS_Cog    ~ 1 + time * group,
          Score_Multi ~ 1 + time * group
        ),
        d = list(~1+time, ~1+time),
        control = list(allow.n.na = TRUE)
      )
    }, error = function(e) return(NULL))
    
    # MMSE models
    # Model 1: Benchmark
    fit_bench <- tryCatch({
      mjoin(
        data = sim_data,
        Y.names = c("MMSE", "AB_PET"), 
        time = "time", 
        long.formulas = list(
          MMSE ~ 1 + time * group,
          AB_PET   ~ 1 + time * group
        ),
        d = list(~1+time, ~1+time),
        control = list(allow.n.na = TRUE) # Allow missing data
      )
    }, error = function(e) return(NULL))

    # Model 2: Multibiomarker
    fit_multi <- tryCatch({
      mjoin(
        data = sim_data,
        Y.names = c("MMSE", "Score_Multi"), 
        time = "time",
        long.formulas = list(
          MMSE    ~ 1 + time * group,
          Score_Multi ~ 1 + time * group
        ),
        d = list(~1+time, ~1+time),
        control = list(allow.n.na = TRUE)
      )
    }, error = function(e) return(NULL))
    
    # --- F. EXTRACT RESULTS ---
    get_trt_est <- function(model, outcome_idx) {
      if(is.null(model)) return(NA)
      # Extract fixed effects
      coefs <- model$coefficients$fixed[[outcome_idx]]
      if("time:group" %in% names(coefs)) return(coefs["time:group"])
      return(NA)
    }
    
    return(data.frame(
      scenario_id = scen$scenario_id,
      replicate   = i,
      
      # ADAS-Cog Results
      bench_ADAS_clin = get_trt_est(fit_bench_adas, 1),
      bench_ADAS_surr = get_trt_est(fit_bench_adas, 2),
      multi_ADAS_clin = get_trt_est(fit_multi_adas, 1),
      multi_ADAS_surr = get_trt_est(fit_multi_adas, 2),
      
      # MMSE Results
      bench_MMSE_clin = get_trt_est(fit_bench_mmse, 1),
      bench_MMSE_surr = get_trt_est(fit_bench_mmse, 2),
      multi_MMSE_clin = get_trt_est(fit_multi_mmse, 1),
      multi_MMSE_surr = get_trt_est(fit_multi_mmse, 2),
      
      N = scen$sample_size,
      Effect = scen$effect_size,
      Corr = scen$corr_strength
    ))

# 6. SAVE FINAL OUTPUT
# Bind the list returned by mclapply into one dataframe
final_results <- bind_rows(results_list)

# Save with Job ID so they don't overwrite each other
output_filename <- paste0("simulation_results_job_", job_id, ".rds")
saveRDS(final_results, output_filename)

cat("Success! Results saved to:", output_filename, "\n")