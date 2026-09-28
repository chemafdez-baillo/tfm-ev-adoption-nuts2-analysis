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

**Regression (M4 — Joint model, n = 236):**

| Variable | β | p-value |
|---|---|---|
| SC1 (socioeconomic composite) | 0.826 | < 0.001 |
| tert_share (tertiarisation) | −1.708 | < 0.001 |
| RC1 (accessible charging supply) | 0.331 | < 0.001 |
| RC2 (charging-to-EV ratio) | −0.145 | < 0.05 |
| own_tax_gdp (property tax exemption) | 53.81 | < 0.01 |
| pol_1yr_gdp (total 1st-year incentive) | −0.835 | 0.075 |

**Four regional profiles (K-means, k = 4, n = 223):**

| Cluster | Label | n | ev_share | Profile |
|---|---|---|---|---|
| C1 | **Pioneers** | 15 | 4.49% | High SC1, very high RC1, low incentives — organic adoption driven by structural conditions. Netherlands, Belgium, Denmark, NW Germany |
| C4 | **Prosperous & Tertiarised** | 27 | 1.43% | High SC1, highest tert_share — adoption below expectations despite high income; advanced-services metros |
| C3 | **European Average** | 148 | 2.11% | Neutral on all dimensions — France, Spain, Italy, Germany, Poland; largest margin for incremental improvement |
| C2 | **Low Development & Compensatory Policy** | 33 | 0.37% | Low SC1, low RC1, highest incentives — financial incentives cannot overcome structural barriers; Eastern and Southern Europe |

**Main conclusions:**
- Socioeconomic factors are the strongest determinant (adj. R² = 0.605 for M1 alone)
- Charging infrastructure has a positive and independent effect on EV adoption
- Financial incentives show a negative coefficient due to endogeneity: countries with lower organic adoption apply more generous incentives precisely to compensate for structural disadvantage

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

## 👤 Author

**José María Fernández-Baillo García-Mascaraque**  
MSc Economic Data Modelling and Analysis — Universidad de Castilla-La Mancha  
📍 Ciudad Real, Spain  
🔗 [LinkedIn](https://www.linkedin.com/in/chema-fernández-baillo-garcía-mascaraque)
