# Aim: Manually assemble ADNI and Literature parameters into the final .rds file

# Library
library(tidyverse)

# 1. Define mean slopes (Fixed Effects / Beta)
# Units: Change per year.

mean_slopes <- c(
  "ADAS_Cog"   = 2.56,   # ADNI
  "MMSE"       = -1.06,  # ADNI
  "AB_PET"     = 2.8,   # Ye et al., 2018
  "tau_PET"    = 0.05,   # Oh et al., 2024
  "vMRI"       = -216.3,   # Sutphen et al., 2018
  "CSF_ptau"   = 0.66, # ADNI
  "CSF_AB42"   = -20.81, # ADNI
  "Plasma_ptau"= 0.14 # ADNI
)

# 2. Define SD of slopes (Random Effects Variability)
# Units: SD of the rate of change.

slope_sds <- c(
  "ADAS_Cog"   = 2.66,   # ADNI
  "MMSE"       = 1.29,   # ADNI
  "AB_PET"     = 5.61,   # Ye et al., 2018
  "tau_PET"    = 0.07,   # Oh et al., 2024
  "vMRI"       = 11.9,   # Sutphen et al., 2018
  "CSF_ptau"   = 1.22,   # ADNI
  "CSF_AB42"   = 26.73,  # ADNI
  "Plasma_ptau"= 5.14    # ADNI  
)

# 3. Define correlations (The D-Matrix Skeleton)
# Define the correlation matrix for the slopes.
# 1.0 on diagonal. Off-diagonals from ADNI or Literature proxies.

# Initialise an 8x8 Identity Matrix
outcomes <- names(mean_slopes)
cor_matrix <- diag(1, nrow = length(outcomes))
dimnames(cor_matrix) <- list(outcomes, outcomes)

#  Fill Cognition <-> Cognition Correlations ---
cor_matrix["ADAS_Cog", "MMSE"] <- -0.902  
cor_matrix["MMSE", "ADAS_Cog"] <- -0.902

# Fill Biomarker <-> Cognition Correlations 
# CSF AB42 (ADNI)
cor_matrix["CSF_AB42", "MMSE"] <- -0.269      
cor_matrix["MMSE", "CSF_AB42"] <- -0.269
cor_matrix["CSF_AB42", "ADAS_Cog"] <- 0.284   
cor_matrix["ADAS_Cog", "CSF_AB42"] <- 0.284
# CSF ptau (ADNI)
cor_matrix["CSF_ptau", "MMSE"] <- 0.027      
cor_matrix["MMSE", "CSF_ptau"] <- 0.027
cor_matrix["CSF_ptau", "ADAS_Cog"] <- -0.039   
cor_matrix["ADAS_Cog", "CSF_ptau"] <- -0.039
# Plasma ptau (ADNI)
cor_matrix["Plasma_ptau", "MMSE"] <- 0.1      
cor_matrix["MMSE", "Plasma_ptau"] <- 0.1
cor_matrix["Plasma_ptau", "ADAS_Cog"] <- -0.093   
cor_matrix["ADAS_Cog", "Plasma_ptau"] <- -0.093
# Tau PET (Singleton)
cor_matrix["tau_PET", "MMSE"] <- -0.56      # From Singleton
cor_matrix["MMSE", "tau_PET"] <- -0.56
cor_matrix["tau_PET", "ADAS_Cog"] <- 0.56   # Derived
cor_matrix["ADAS_Cog", "tau_PET"] <- 0.56
# vMRI (Singleton)
cor_matrix["vMRI", "MMSE"] <- 0.55          # Thickness goes DOWN, MMSE goes DOWN (Positive corr)
cor_matrix["MMSE", "vMRI"] <- 0.55         
cor_matrix["vMRI", "ADAS_Cog"] <- -0.55     # Thickness DOWN, ADAS UP (Negative corr)
cor_matrix["ADAS_Cog", "vMRI"] <- -0.55
# Aβ-PET (Budd Haeberlein) ["Weak" Anchor]
cor_matrix["AB_PET", "MMSE"] <- -0.24
cor_matrix["MMSE", "AB_PET"] <- -0.24
cor_matrix["AB_PET", "ADAS_Cog"] <- 0.20
cor_matrix["ADAS_Cog", "AB_PET"] <- 0.20


# Fill Biomarker <-> Biomarker Correlations
cor_matrix["CSF_AB42", "CSF_ptau"] <- 0.004 # ADNI      
cor_matrix["CSF_ptau", "CSF_AB42"] <- 0.004
cor_matrix["CSF_AB42", "AB_PET"] <- -0.4 # SMALL PROXY      
cor_matrix["AB_PET", "CSF_AB42"] <- -0.4
cor_matrix["CSF_AB42", "tau_PET"] <- -0.45 # Literature     
cor_matrix["tau_PET", "CSF_AB42"] <- -0.45
cor_matrix["CSF_AB42", "vMRI"] <- 0.1 # SMALL PROXY     
cor_matrix["vMRI", "CSF_AB42"] <- 0.1
cor_matrix["CSF_AB42", "Plasma_ptau"] <- -0.031 # ADNI     
cor_matrix["Plasma_ptau", "CSF_AB42"] <- -0.031
cor_matrix["CSF_ptau", "AB_PET"] <- 0.1 # SMALL PROXY     
cor_matrix["AB_PET", "CSF_ptau"] <- 0.1
cor_matrix["CSF_ptau", "tau_PET"] <- 0.5 # PROXY     
cor_matrix["tau_PET", "CSF_ptau"] <- 0.5
cor_matrix["CSF_ptau", "vMRI"] <- -0.37 # PROXY     
cor_matrix["vMRI", "CSF_ptau"] <- -0.37
cor_matrix["CSF_ptau", "Plasma_ptau"] <- 0.047 # ADNI     
cor_matrix["Plasma_ptau", "CSF_ptau"] <- 0.047
cor_matrix["AB_PET", "Plasma_ptau"] <- 0.36 # Literature     
cor_matrix["Plasma_ptau", "AB_PET"] <- 0.36
cor_matrix["tau_PET", "Plasma_ptau"] <- 0.4 # Literature     
cor_matrix["Plasma_ptau", "tau_PET"] <- 0.4
cor_matrix["vMRI", "Plasma_ptau"] <- -0.106 # Literature     
cor_matrix["Plasma_ptau", "vMRI"] <- -0.106

# 4. Construct final objects

# Convert Correlation Matrix to Covariance Matrix (D-Matrix)
# Formula: Cov(x,y) = Cor(x,y) * SD(x) * SD(y) (Helwig et al., 2017)
sd_matrix <- diag(slope_sds)
D_matrix <- sd_matrix %*% cor_matrix %*% sd_matrix
dimnames(D_matrix) <- list(outcomes, outcomes)

# Define Intercepts (ADNI baseline)
intercepts <- c(
  "ADAS_Cog" = 18.5, "MMSE" = 26, "AB_PET" = 79.8, # Centiloids, 
  "tau_PET" = 1.5, "vMRI" = 7375.4, "CSF_ptau" = 31.7,
  "CSF_AB42" = 793.9, "Plasma_ptau" = 22.1
)

# Define Residual Errors (Sigma)
residuals <- c(
  "ADAS_Cog" = 3.86, "MMSE" = 1.75 , "AB_PET" = 10, 
  "tau_PET" = 0.15 , "vMRI" = 100, "CSF_ptau" = 2.89,
  "CSF_AB42" = 124.31 , "Plasma_ptau" = 8.78
)

# 5. Save
simulation_parameters <- list(
  beta = mean_slopes,      # Vector of mean slopes
  intercepts = intercepts, # Vector of mean baselines
  D_matrix = D_matrix,     # Covariance matrix of slopes
  sigma = residuals,       # Vector of residual SDs
  cor_matrix_base = cor_matrix # Save the base correlations for reference
)

# Save to .rds
saveRDS(simulation_parameters, "simulation_parameters.rds")

cat("--- SUCCESS ---\n")
cat("simulation_parameters.rds has been created.\n")
cat("Upload this file to Eddie.\n")
