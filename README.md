# MIIPW: GEE with Missing Data Using Mean Score, IPW, and Multiple Imputation
<img src="man/figures/MIIPWlogo.png" align="right" alt="MIIPW logo" width="200">

<!-- badges: start -->

[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![CRAN status](https://www.r-pkg.org/badges/version/MIIPW)](https://CRAN.R-project.org/package=MIIPW)
[![R build status](https://github.com/yourusername/MIIPW/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/MIIPW/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

<br/>

> **MIIPW** is an R package for fitting Generalized Estimating Equations (GEE) when data is missing at random (MAR), using a combination of **Mean Score (MS)**, **Inverse Probability Weighting (IPW)**, and **Multiple Imputation (MI)** techniques.

---

## ✨ Features

The MIIPW package supports robust marginal modeling with missing outcome and covariate data by integrating:

- 📌 **Mean Score Approach** for consistent estimation under missingness
- 🧮 **IPW Estimation** using modeled response or missingness probabilities
- ♻️ **Multiple Imputation (MI)** integration for general missing data structures
- 📈 **Flexible model specification** with various correlation structures
- 📦 Easy interface using `gee()` style model formulas

---

## 📦 Installation

### Via `pak` (Recommended)

```r
install.packages("pak")
pak::pak("yourusername/MIIPW")
