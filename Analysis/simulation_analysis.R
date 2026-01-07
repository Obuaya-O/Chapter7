# Library
library(tidyverse)

# 1. Load the 4 Raw Datasets
res_w_adas <- readRDS("FINAL_weighted_ADAS.rds") %>% mutate(Model = "Weighted_ADAS")
res_w_mmse <- readRDS("FINAL_weighted_MMSE.rds") %>% mutate(Model = "Weighted_MMSE")
res_u_adas <- readRDS("FINAL_unweighted_ADAS.rds") %>% mutate(Model = "Unweighted_ADAS")
res_u_mmse <- readRDS("FINAL_unweighted_MMSE.rds") %>% mutate(Model = "Unweighted_MMSE")

# Sanity check
colnames(res_u_adas)

# Combine into one master list
all_data <- bind_rows(res_w_adas, res_w_mmse, res_u_adas, res_u_mmse)

# 2. Define True Effects (for Bias calculation)
# Must match setup_params.R
get_true_effect <- function(eff_size, outcome) {
  # These slopes come from simulation constants
  if (eff_size == "Null") return(0)
  
  # Multiplier was -0.25 (Medium) or -0.50 (Large)
  mod <- ifelse(eff_size == "Medium", -0.25, -0.50)
  
  if (outcome == "ADAS") return(2.56 * mod) # SLOPE_ADAS * Modifier
  if (outcome == "MMSE") return(-1.06 * mod) # SLOPE_MMSE * Modifier
  return(NA)
}

# 3. Calculate Performance Metrics (Bias, MSE, ESE)
performance_summary <- all_data %>%
  group_by(Model, N, Corr, Miss, effect_size) %>% # Group by Scenario
  summarise(
    # True Effect for this specific group
    True_Clin = get_true_effect(unique(effect_size), ifelse(grepl("ADAS", unique(Model)), "ADAS", "MMSE")),
    
    # --- Bias (Mean Estimate - True Effect) ---
    Bias_Clin = mean(est_bench_clin, na.rm=TRUE) - True_Clin,
    Bias_Surr = mean(est_multi_surr, na.rm=TRUE) - 0, # Assuming Null surrogate effect if not defined, or calculate true surrogate slope
    
    # --- ESE (Empirical Standard Error = SD of Estimates) ---
    ESE_Clin = sd(est_bench_clin, na.rm=TRUE),
    ESE_Surr = sd(est_multi_surr, na.rm=TRUE),
    
    # --- MSE (Bias^2 + Variance) ---
    MSE_Clin = Bias_Clin^2 + var(est_bench_clin, na.rm=TRUE),
    
    .groups = "drop"
  )

print("Performance Metrics Calculated.")

# R2 and STE
# Function to calculate R2_trial and STE from a set of 1000 trials
# Function to calculate R2 and STE for BOTH methods
calc_surrogacy_comparison <- function(df) {
  
  # 1. Benchmark Model (Amyloid -> Clinical)
  fit_bench <- lm(est_bench_clin ~ est_bench_surr, data = df)
  r2_bench  <- summary(fit_bench)$r.squared
  
  # 2. Multibiomarker Model (Panel -> Clinical)
  fit_multi <- lm(est_bench_clin ~ est_multi_surr, data = df)
  r2_multi  <- summary(fit_multi)$r.squared
  
  # Return one row with both values
  return(data.frame(
    R2_Bench = r2_bench,
    R2_Multi = r2_multi
  ))
}

# Apply to every scenario
surrogacy_summary <- all_data %>%
  group_by(Model, N, Corr, Miss, effect_size) %>%
  do(calc_surrogacy_comparison(.)) %>%
  ungroup()

print("Surrogacy (R2) Calculated for both Benchmark and Multi.")

# FIGURE 1: Surrogacy validation
# Objective: Compare R2_trial of Benchmark vs. Unweighted vs. Weighted Multibiomarker

# Create 3 separate dataframes (one for each bar type) and stack them.

# A) Benchmark Data (We only need to take this from one model to avoid duplicates)
df_bench <- surrogacy_summary %>%
  filter(effect_size != "Null", grepl("Weighted", Model)) %>% # Filter to just Weighted rows to grab unique benchmarks
  select(N, Corr, Miss, R2 = R2_Bench) %>%
  mutate(Type = "Benchmark (Amyloid)")

# B) Weighted Panel Data
df_weighted <- surrogacy_summary %>%
  filter(effect_size != "Null", grepl("Weighted", Model)) %>%
  select(N, Corr, Miss, R2 = R2_Multi) %>%
  mutate(Type = "Weighted Panel")

# C) Unweighted Panel Data
df_unweighted <- surrogacy_summary %>%
  filter(effect_size != "Null", grepl("Unweighted", Model)) %>%
  select(N, Corr, Miss, R2 = R2_Multi) %>%
  mutate(Type = "Unweighted Panel")

# Combine them all
plot_ready <- bind_rows(df_bench, df_unweighted, df_weighted)

# Set Order: Benchmark -> Unweighted -> Weighted
plot_ready$Type <- factor(plot_ready$Type, 
                          levels = c("Benchmark (Amyloid)", "Unweighted Panel", "Weighted Panel"))

# 2. Plot
p1 <- ggplot(plot_ready, aes(x = Corr, y = R2, fill = Type)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.7) +
  # Split by Missing Data (Rows) and Sample Size (Cols)
  facet_grid(Miss ~ N, labeller = label_both) + 
  theme_bw(base_size = 11) +
  labs(
    title = "Surrogacy Strength Validation",
    subtitle = "Comparison: Benchmark vs. Unweighted vs. Weighted Multibiomarker",
    y = expression(R[trial]^2 ~ "(Correlation Strength)"),
    x = "Correlation Structure",
    fill = "Method"
  ) +
  # Define your 3 colors here
  scale_fill_manual(values = c(
    "Benchmark (Amyloid)" = "gray70", 
    "Unweighted Panel"    = "#7fcdbb",  # Teal/Green 
    "Weighted Panel"      = "#2c7fb8"   # Strong Blue
  )) +
  theme(legend.position = "bottom")

print(p1)

# Get exact numbers because graphs are just for looks
r2_table <- plot_ready %>%
  group_by(Type, N, Corr, Miss) %>%
  summarise(
    Mean_R2 = mean(R2, na.rm = TRUE),
    SD_R2   = sd(R2, na.rm = TRUE), # Optional: Standard Deviation to see variability
    .groups = "drop"
  ) %>%
  # 2. Reshape for readability (Methods side-by-side)
  pivot_wider(
    names_from = Type, 
    values_from = c(Mean_R2, SD_R2)
  ) %>%
  # 3. Clean up columns (Select means mostly, maybe drop SDs for the main view)
  select(
    Sample_Size = N,
    Correlation = Corr,
    Missingness = Miss,
    # Select the Mean columns for the 3 methods
    Bench_Amyloid = `Mean_R2_Benchmark (Amyloid)`,
    Unweighted    = `Mean_R2_Unweighted Panel`,
    Weighted      = `Mean_R2_Weighted Panel`
  ) %>%
  # 4. Calculate the "Win" margin (How much better is Weighted?)
  mutate(
    Gain_vs_Bench = Weighted - Bench_Amyloid,
    Gain_vs_Unweighted = Weighted - Unweighted
  ) %>%
  arrange(Sample_Size, Correlation, Missingness)

# 5. View the Table
print(r2_table)

# 6. Save to CSV for your thesis (Excel)
write.csv(r2_table, "Table_R2_Comparison_Exact_Numbers.csv", row.names = FALSE)

# ESE
# 1. Prepare ESE Data for all 3 Methods
# Note: We Standardise Benchmark by 20 (Simulated Pop SD) to make it comparable to Z-scores

# A) Benchmark ESE (Standardized)
df_ese_bench <- all_data %>%
  filter(effect_size != "Null", grepl("Weighted", Model)) %>% # Take from one model to avoid duplicates
  group_by(N, Corr, Miss) %>%
  summarise(ESE = sd(est_bench_surr, na.rm=TRUE) / 20, .groups="drop") %>% # /20 for Scaling
  mutate(Method = "Benchmark (Amyloid)")

# B) Weighted Panel ESE
df_ese_weight <- all_data %>%
  filter(effect_size != "Null", grepl("Weighted", Model)) %>%
  group_by(N, Corr, Miss) %>%
  summarise(ESE = sd(est_multi_surr, na.rm=TRUE), .groups="drop") %>%
  mutate(Method = "Weighted Panel")

# C) Unweighted Panel ESE
df_ese_unweight <- all_data %>%
  filter(effect_size != "Null", grepl("Unweighted", Model)) %>%
  group_by(N, Corr, Miss) %>%
  summarise(ESE = sd(est_multi_surr, na.rm=TRUE), .groups="drop") %>%
  mutate(Method = "Unweighted Panel")

# Combine
plot_ese_ready <- bind_rows(df_ese_bench, df_ese_weight, df_ese_unweight)

# Set Factor Order
plot_ese_ready$Method <- factor(plot_ese_ready$Method, 
                                levels = c("Benchmark (Amyloid)", "Unweighted Panel", "Weighted Panel"))

# 2. Plot
p2 <- ggplot(plot_ese_ready, aes(x = as.factor(N), y = ESE, group = Method, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_grid(Miss ~ Corr, labeller = label_both) +
  theme_bw(base_size = 11) +
  labs(
    title = "Precision of Surrogate Estimates (ESE)",
    subtitle = "Lower ESE = Higher Precision (Benchmark standardised for comparison)",
    y = "Empirical Standard Error (SD Units)",
    x = "Sample Size (N)",
    color = "Method"
  ) +
  scale_color_manual(values = c(
    "Benchmark (Amyloid)" = "gray60", 
    "Unweighted Panel"    = "#7fcdbb", 
    "Weighted Panel"      = "#2c7fb8"
  )) +
  theme(legend.position = "bottom")

print(p2)

# Get numbers
# 1. Calculate Signal-to-Noise Ratio (SNR)
# SNR = |Mean Estimate| / SD(Estimate)
# This cancels out the units (Centiloids vs Z-scores), allowing direct comparison.

snr_table <- all_data %>%
  filter(effect_size != "Null") %>%
  group_by(N, Corr, Miss, Model) %>%
  summarise(
    # BENCHMARK (Amyloid)
    Signal_Bench = mean(est_bench_surr, na.rm=TRUE),
    Noise_Bench  = sd(est_bench_surr, na.rm=TRUE),
    SNR_Bench    = abs(Signal_Bench) / Noise_Bench,
    
    # PANEL (Multibiomarker)
    Signal_Panel = mean(est_multi_surr, na.rm=TRUE),
    Noise_Panel  = sd(est_multi_surr, na.rm=TRUE),
    SNR_Panel    = abs(Signal_Panel) / Noise_Panel,
    
    .groups = "drop"
  ) %>%
  # 2. Label Methods
  mutate(Method = ifelse(grepl("Weighted", Model), "Weighted", "Unweighted")) %>%
  
  # 3. Clean and Pivot for Comparison
  select(N, Corr, Miss, Method, SNR_Bench, SNR_Panel) %>%
  group_by(N, Corr, Miss) %>%
  summarise(
    # Benchmark is identical for both weighted/unweighted runs, just average them
    Benchmark_SNR = mean(SNR_Bench),
    Unweighted_SNR = mean(SNR_Panel[Method == "Unweighted"]),
    Weighted_SNR   = mean(SNR_Panel[Method == "Weighted"]),
    .groups = "drop"
  ) %>%
  # 4. Calculate the "Win" %
  mutate(
    # Positive % means Weighted is BETTER than Benchmark
    Pct_Improvement_vs_Bench = (Weighted_SNR - Benchmark_SNR) / Benchmark_SNR * 100,
    
    # Positive means Weighted is BETTER than Unweighted
    Gain_vs_Unweighted = Weighted_SNR - Unweighted_SNR
  ) %>%
  arrange(N, Corr, Miss)

# 5. View and Save
print(head(snr_table))
write.csv(snr_table, "Table_SNR_Power_Comparison.csv", row.names = FALSE)

# 1. Use the SNR table you just created
power_table <- snr_table %>%
  mutate(
    # Formula: pnorm( SNR - Z_score_cutoff )
    # 1.96 is the Z-cutoff for p < 0.05 (two-sided)
    
    Power_Bench_Pct    = pnorm(Benchmark_SNR - 1.96) * 100,
    Power_Unweighted_Pct = pnorm(Unweighted_SNR - 1.96) * 100,
    Power_Weighted_Pct   = pnorm(Weighted_SNR - 1.96) * 100
  ) %>%
  select(N, Corr, Miss, 
         SNR_Bench = Benchmark_SNR, Power_Bench_Pct,
         SNR_Weight = Weighted_SNR, Power_Weighted_Pct)

# 2. View the "Implied Power"
print(head(power_table))
write.csv(power_table, "Table_Theoretical_Power.csv", row.names = FALSE)

# BIAS
# 1. Define True Clinical Effects function (Simulation Constants)
get_true_clin <- function(eff, outcome) {
  if (eff == "Null") return(0)
  mod <- ifelse(eff == "Medium", -0.25, -0.50)
  if (outcome == "ADAS") return(2.56 * mod) 
  if (outcome == "MMSE") return(-1.06 * mod)
  return(NA)
}

# 2. Prepare Bias Data for ALL 3 Methods
df_bias_all <- all_data %>%
  filter(effect_size != "Null") %>%
  mutate(
    Outcome = ifelse(grepl("ADAS", Model), "ADAS", "MMSE"),
    Method = case_when(
      grepl("Weighted", Model) ~ "Weighted Panel",
      grepl("Unweighted", Model) ~ "Unweighted Panel",
      TRUE ~ "Benchmark (Amyloid)" # Fallback, though usually derive benchmark from weighted file
    )
  ) %>%
  # IMPORTANT: Make sure to treat the Benchmark properly.
  # Since the 'Model' column distinguishes files, group by Model 
  # and then rename them to the 3 Types.
  group_by(N, Corr, Miss, effect_size, Outcome, Model) %>%
  summarise(
    Mean_Est = mean(est_bench_clin, na.rm = TRUE),
    True_Val = get_true_clin(unique(effect_size), unique(Outcome)),
    Bias = Mean_Est - True_Val,
    .groups = "drop"
  ) %>%
  mutate(
    Method_Label = case_when(
      grepl("Weighted", Model) ~ "Weighted Panel",
      grepl("Unweighted", Model) ~ "Unweighted Panel",
      # Note: Benchmark data is inside the weighted/unweighted files, 
      # clinical estimates are identical. Just plot the two panels 
      # as they represent the simulation runs.
      TRUE ~ "Benchmark"
    )
  )

# 3. Plotting
# Facet by Outcome (ADAS vs MMSE) to keep scales correct
p3 <- ggplot(df_bias_all, aes(x = as.factor(N), y = Bias, fill = Method_Label)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Facet: Rows = Missing Data, Cols = Effect Size
  facet_grid(Miss ~ effect_size, labeller = label_both) +
  theme_bw(base_size = 11) +
  labs(
    title = "Estimation Bias (Clinical Outcome)",
    subtitle = "Difference between Estimated and True Treatment Effect (Target = 0)",
    y = "Bias (Estimate - Truth)",
    x = "Sample Size (N)",
    fill = "Simulation Run"
  ) +
  scale_fill_manual(values = c("Weighted Panel" = "#2c7fb8", "Unweighted Panel" = "#7fcdbb")) +
  coord_flip()

print(p3)

# Exact numbers
# 1. Generate the Table
bias_table <- df_bias_all %>%
  # KEEP effect_size here so rows remain unique
  select(N, Corr, Miss, Outcome, effect_size, Method_Label, Bias) %>%
  
  # Pivot to see methods side-by-side
  pivot_wider(names_from = Method_Label, values_from = Bias) %>%
  
  # Sort for clean reading
  arrange(Outcome, N, Corr, Miss, effect_size)

# 2. Check structure (Should show 'num' columns, not 'list')
str(bias_table)

# 2. View and Save
print(head(bias_table))
write.csv(bias_table, "Table_Bias_Clinical_Exact_Numbers.csv", row.names = FALSE)
