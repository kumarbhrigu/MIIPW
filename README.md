# MIIPW

**MIIPW** (Mean Score and Inverse Probability Weighted Methods with Multiple Imputation) is an R package for analyzing longitudinal or repeated measures data with missing responses or covariates. It implements advanced methods such as inverse probability weighting (IPW), mean score estimation, and multiple imputation to provide robust parameter estimates using GEE (Generalized Estimating Equations).

## Key Features

- Implements **Augmented IPW (AIPW)** and **Simple IPW (SIPW)** for longitudinal GEE models.
- Provides **Mean Score Estimation** with multiple imputation to address missing data in covariates or response.
- Includes functions for fitting **marginal models** using semiparametric methods.
- Supports various correlation structures: `"independence"`, `"exchangeable"`, `"AR-1"`, and `"unstructured"`.
- Functions for **model selection** using QIC (Quasi-likelihood under the Independence Model Criterion).
- Includes internal tools for parameter estimation using **Fisher Scoring**.

## Installation

```r
# Install from GitHub (requires remotes or devtools)
# remotes::install_github("kumarbhrigu/MIIPW")
