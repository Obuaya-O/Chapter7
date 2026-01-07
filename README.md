# Multibiomarker surrogate validation simulation

## Overview
This repository contains the R code and simulation scripts for **Chapter 7: Developing and evaluating a model to identify and prioritise candidate surrogate outcomes for Alzheimer's disease clinical trials**.

The purpose of this study was to evaluate whether a Multibiomarker Panel (combining amyloid PET, tau PET, vMRI and fluid biomarkers) provides superior predictive validity for cognitive decline compared to a single-biomarker benchmark (amyloid PET alone). Using a Monte Carlo simulation framework, we generated synthetic clinical trial data to assess surrogate performance under varying conditions of sample size, effect size and missing data mechanisms.

## Key Findings
* **Surrogacy:** The weighted multibiomarker panel achieved a trial-level $R^2 \approx 0.25$, significantly outperforming the single-biomarker benchmark ($R^2 < 0.05$).
* **Efficiency:** Despite higher measurement noise, the composite score improved the Signal-to-Noise Ratio (SNR), increasing statistical power from ~45% to ~77% in small trials ($N=260$).
* **Robustness:** The model was robust to Missing at Random (MAR) assumptions but displayed significant performance degradation (~27% drop in $R^2$) under informative censoring (MNAR).

## Repository Structure

### 1. Data Generation 
Scripts to generate synthetic subject-level data based on parameters derived from the Alzheimer's Disease Neuroimaging Initiative and Phase III clinical trials.

### 2. Analysis 
Scripts to fit Linear Mixed-Effects Models (`lme4`) and calculate trial-level surrogacy ($R^2_{trial}$).

### 3. Results 
Full result files

## Dependencies
The simulation relies on the following R packages:
* **`lme4`**: For fitting linear mixed-effects models.
* **`MASS`**: For generating multivariate normal distributions (`mvrnorm`).
* **`tidyverse`**: For data manipulation and visualisation.

## References and Inspiration

This simulation framework was informed by best practices and methodologies from the following key sources:

* **Simulation Design:**
    * *Morris, T. P., White, I. R., & Crowther, M. J. (2019).* Using simulation studies to evaluate statistical methods. *Statistics in Medicine*. [Link](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8086)
* **Statistical Modelling (LMM):**
    * *Bates, D., et al. (2015).* Fitting Linear Mixed-Effects Models Using lme4. *Journal of Statistical Software*. [Link](https://www.jstatsoft.org/article/view/v067i01)
    * *Murphy, J. I., Weaver, N. E., & Hendricks, A. E. (2022).* Accessible analysis of longitudinal data with linear mixed effects models. *Disease Models & Mechanisms*. [Link](https://journals.biologists.com/dmm/article/15/5/dmm048025/275308/Accessible-analysis-of-longitudinal-data-with)

* **Multivariate Data Generation:**
    * *Moore, R., et al. (2019).* A linear mixed-model approach to study multivariate gene-environment interactions. *Nature Genetics*. (Methodological reference for handling multivariate correlated components). [Link](https://www.nature.com/articles/s41588-018-0271-0)
