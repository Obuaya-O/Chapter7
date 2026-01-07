#!/usr/bin/env Rscript
# PURPOSE: Executes the Monte Carlo simulation to evaluate the Unweighted Biomarker Panel using ADAS-Cog as the primary clinical endpoint

lib_path <- "/Onyeka_Rlib"
.libPaths(lib_path)

# Libraries
library(dplyr); library(tidyr); library(MASS); library(Matrix); library(lme4) #tidyverse was crashing for me so I loaded only what I needed

get_trt_est <- function(model) {
  if(is.null(model)) return(NA)
  tryCatch({
    coefs <- summary(model)$coefficients
    if("time:group" %in% rownames(coefs)) return(coefs["time:group", "Estimate"])
    return(NA)
  }, error = function(e) NA)
}

args <- commandArgs(trailingOnly = TRUE)
job_id <- if (length(args) == 0) 1 else as.integer(args[1])
scenarios <- readRDS("simulation_parameters.rds")
if(job_id > nrow(scenarios)) quit(save = "no")
scen <- scenarios[job_id, ]
print(paste("Unweighted ADAS Scen:", job_id, "| N:", scen$sample_size))

# Constants
SLOPE_ADAS <- 2.56
SLOPE_MMSE <- -1.06
S_AB   <- 2.8
S_TAU  <- -216.3
S_VMRI <- 0.05
S_CSFP <- 0.66
S_CSFA <- -20.81
S_PLAS <- 0.14

TRT_MOD <- switch(as.character(scen$effect_size), "Null"=0, "Medium"=-0.25, "Large"=-0.50)

chunk_id <- scen$chunk_id # Had to chunk the results due to crashing issues on the cluster
reps_per_chunk <- 50 # (1000 total / 20 chunks)
results <- data.frame()

for (i in 1:reps_per_chunk) {
  # 1. COVARIANCE
  D <- matrix(0, 7, 7)
  diag(D) <- c(1.0, 2.0, 0.56, 0.55, 0.039, 0.284, 0.093)
  
  if (scen$corr_level == "Real") rhos <- c(0.20, 0.56, 0.55, -0.027, 0.269, -0.100)
  else if (scen$corr_level == "Moderate") rhos <- rep(0.4, 6)
  else rhos <- rep(0.6, 6)
  
  D[1,2] <- D[2,1] <- rhos[1] * sqrt(D[1,1] * D[2,2])
  D[1,3] <- D[3,1] <- rhos[2] * sqrt(D[1,1] * D[3,3])
  D[1,4] <- D[4,1] <- rhos[3] * sqrt(D[1,1] * D[4,4])
  D[1,5] <- D[5,1] <- rhos[4] * sqrt(D[1,1] * D[5,5])
  D[1,6] <- D[6,1] <- rhos[5] * sqrt(D[1,1] * D[6,6])
  D[1,7] <- D[7,1] <- rhos[6] * sqrt(D[1,1] * D[7,7])

  ev <- eigen(D)$values; if(any(ev<=1e-6)) D <- as.matrix(nearPD(D)$mat)
  
  # START GENERATION BLOCK 
  n_gen <- scen$sample_size * 2
  b_i_gen <- mvrnorm(n_gen, mu=rep(0,7), Sigma=D)
  colnames(b_i_gen) <- c("b_a","b_ab","b_tau","b_vmri","b_cp","b_ca","b_pl")
  
  age_gen <- rnorm(n_gen, mean=73, sd=7.5)
  
  subjects_pool <- data.frame(sub_id_gen = 1:n_gen, Age = age_gen, group = sample(0:1, n_gen, replace=TRUE))
  subjects_pool <- cbind(subjects_pool, b_i_gen)
  
  subjects_pool$ADAS_BL   <- 15  + rnorm(n_gen, 0, 8) 
  subjects_pool$MMSE_BL   <- 26  + rnorm(n_gen, 0, 2)
  subjects_pool$AB_PET_BL <- 60  + rnorm(n_gen, 0, 20) 
  
  valid_pool <- subjects_pool %>%
    filter(Age >= 50 & Age <= 90) %>%
    filter(MMSE_BL >= 16 & MMSE_BL <= 28) %>%
    filter(ADAS_BL >= 3.9 & ADAS_BL <= 62.9) %>%
    filter(AB_PET_BL >= 24)
  
  if(nrow(valid_pool) < scen$sample_size) stop("Increase oversampling factor.")
  
  selected_subs <- valid_pool[1:scen$sample_size, ]
  selected_subs$sub_id <- 1:scen$sample_size 
  
  sim_d <- expand.grid(sub_id=1:scen$sample_size, time=seq(0,2,0.5)) %>% left_join(selected_subs, by="sub_id")
  grp <- sim_d$group; t <- sim_d$time
  
  sim_d$ADAS <- sim_d$ADAS_BL + (SLOPE_ADAS*(1+TRT_MOD*grp) + sim_d$b_a) * t + rnorm(nrow(sim_d), 0, 1)
  sim_d$MMSE <- sim_d$MMSE_BL + (SLOPE_MMSE*(1+TRT_MOD*grp) + sim_d$b_a) * t + rnorm(nrow(sim_d), 0, 1)
  sim_d$AB_PET <- sim_d$AB_PET_BL + (S_AB*(1+TRT_MOD*grp) + sim_d$b_ab) * t + rnorm(nrow(sim_d), 0, 5)
  sim_d$Tau    <- (200 + rnorm(nrow(sim_d),0,40)) + (S_TAU*(1+TRT_MOD*grp)  + sim_d$b_tau) * t + rnorm(nrow(sim_d),0,10) 
  sim_d$vMRI   <- (1.2 + rnorm(nrow(sim_d),0,0.2)) + (S_VMRI*(1+TRT_MOD*grp) + sim_d$b_vmri)* t + rnorm(nrow(sim_d),0,0.1)
  sim_d$CSFP   <- (1.2 + rnorm(nrow(sim_d),0,0.2)) + (S_CSFP*(1+TRT_MOD*grp) + sim_d$b_cp)  * t + rnorm(nrow(sim_d),0,0.1)
  sim_d$CSFA   <- (200 + rnorm(nrow(sim_d),0,40)) + (S_CSFA*(1+TRT_MOD*grp)  + sim_d$b_ca)  * t + rnorm(nrow(sim_d),0,10)
  sim_d$Plas   <- (1.2 + rnorm(nrow(sim_d),0,0.2)) + (S_PLAS*(1+TRT_MOD*grp) + sim_d$b_pl)  * t + rnorm(nrow(sim_d),0,0.1)

  if (scen$missing_type == "MCAR") sim_d <- sim_d %>% filter(time == 0 | runif(n()) > 0.15)
  else if (scen$missing_type == "MAR") sim_d <- sim_d %>% filter(time == 0 | ifelse(Age > 75, runif(n()) > 0.20, runif(n()) > 0.10))
  else if (scen$missing_type == "MNAR") {
    fast <- which(selected_subs$b_a > quantile(selected_subs$b_a, 0.8))
    sim_d <- sim_d %>% filter(!(sub_id %in% fast & time == 2.0))
  }
  
  # PANEL CALCULATION (Unweighted) 
  std <- function(x) (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  w <- rep(1/6, 6)
  
  sim_d$W_Panel <- (w[1]*std(sim_d$AB_PET) - w[2]*std(sim_d$Tau) + w[3]*std(sim_d$vMRI) + 
                  w[4]*std(sim_d$CSFP)   - w[5]*std(sim_d$CSFA)+ w[6]*std(sim_d$Plas))

  # MODEL FITS (ADAS Target)
  fit_clin  <- tryCatch(lmer(ADAS ~ time * group + (1 | sub_id), data = sim_d), error=function(e) NULL)
  fit_ab    <- tryCatch(lmer(AB_PET ~ time * group + (1 | sub_id), data = sim_d), error=function(e) NULL)
  fit_multi <- tryCatch(lmer(W_Panel ~ time * group + (1 | sub_id), data = sim_d), error=function(e) NULL)

  res_row <- data.frame(
    iter = i,
    est_bench_clin = get_trt_est(fit_clin),
    est_multi_clin = get_trt_est(fit_clin),
    est_bench_surr = get_trt_est(fit_ab),
    est_multi_surr = get_trt_est(fit_multi),
    N = scen$sample_size, Corr = scen$corr_level, Miss = scen$missing_type
  )
  results <- rbind(results, res_row)
}

dir.create("results_unweighted_ADAS", showWarnings=F)
saveRDS(results, paste0("results_unweighted_ADAS/job_", job_id, "_chunk_", chunk_id, ".rds"))