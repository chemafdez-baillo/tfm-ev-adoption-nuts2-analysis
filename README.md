# 🔋 Electric Vehicle Adoption Across European NUTS-2 Regions
### A multivariate econometric analysis of regional EV adoption patterns in the EU

---

## 📋 Overview

This project analyses the determinants of electric vehicle (EV) adoption across 
**236 NUTS-2 regions from 26 EU countries**, combining dimensionality reduction, 
multivariate regression and clustering techniques to identify distinct regional 
adoption profiles.

Developed as a Master's thesis for the **MUMADE programme** 
(Economic Data Modelling and Analysis) at the Universidad de Castilla-La Mancha, Spain.

---

## 🎯 Research Objective

To identify which socioeconomic, infrastructural and geographic factors 
best explain the heterogeneity in EV adoption rates across European regions, 
and to classify regions into interpretable profiles for policy and business targeting.

---

## 🔬 Methodology

The analysis follows a four-stage pipeline:

Raw data (XLSX)
│
▼

EXPLORATORY ANALYSIS
· Correlation matrix (corrplot, ggcorrplot)
· Multicollinearity diagnostics
· KMO & Bartlett's test of sphericity (psych)
│
▼
DIMENSIONALITY REDUCTION
· Principal Component Analysis (PCA)
· Varimax rotation for interpretability (factoextra)
│
▼
MULTIVARIATE REGRESSION
· OLS with HC3 robust standard errors (sandwich, lmtest)
· Model diagnostics: VIF, Breusch-Pagan, residual analysis (car, broom)
· Adjusted R² = 0.676
│
▼
CLUSTERING
· K-means clustering → 4 regional profiles (cluster, factoextra)
· Silhouette analysis for optimal k selection


---

## 📊 Key Findings

- The OLS model with robust errors explains **67.6% of the variance** 
  in regional EV adoption rates (Adjusted R² = 0.676).
- Four distinct regional profiles were identified:
  - **Profile 1** — High adoption: wealthy urban regions with dense charging 
    infrastructure (Nordic countries, NL, AT)
  - **Profile 2** — Growing markets: mid-income regions with strong policy 
    incentives (FR, DE, ES metropolitan areas)
  - **Profile 3** — Lagging regions: lower-income areas with sparse 
    infrastructure (Eastern EU)
  - **Profile 4** — Rural/peripheral: low density, limited grid access, 
    minimal adoption regardless of income
- Charging infrastructure density and regional GDP per capita emerged as 
  the strongest predictors across all model specifications.

---

## 🛠️ Tech Stack

| Category | Packages |
|---|---|
| Data import | `readxl` |
| Data wrangling | `dplyr`, `tidyr`, `scales` |
| Correlation & EDA | `corrplot`, `ggcorrplot`, `PerformanceAnalytics`, `psych` |
| PCA & clustering | `factoextra`, `cluster`, `dbscan` |
| Visualisation | `ggplot2`, `patchwork` |
| Regression & diagnostics | `lmtest`, `sandwich`, `car`, `broom` |
| Reporting | `kableExtra` |

---

## 📁 Repository Structure

📦 tfm-ev-adoption-europe
┣ 📂 data
┃ ┗ 📄 nuts2_ev_dataset.xlsx # Raw dataset (236 regions × variables)
┣ 📂 scripts
┃ ┣ 📄 01_eda_correlations.R # Exploratory analysis & correlation matrices
┃ ┣ 📄 02_pca_analysis.R # PCA with varimax rotation
┃ ┣ 📄 03_regression_ols.R # OLS regression with robust errors
┃ ┗ 📄 04_clustering_kmeans.R # K-means clustering & profiling
┣ 📂 outputs
┃ ┣ 📂 figures # Exported plots (PNG)
┃ ┗ 📂 tables # Exported regression & cluster tables
┗ 📄 README.md


---

## 👤 Author

**José María Fernández-Baillo García-Mascaraque**  
MSc Economic Data Modelling and Analysis — Universidad de Castilla-La Mancha  
📍 Ciudad Real, Spain  
🔗 [LinkedIn](https://www.linkedin.com/in/chema-fernández-baillo-garcía-mascaraque)
