# HIV Immunological Evaluation

This repository contains the scripts and datasets used in the following studies:

1. **“sVCAM-1 and hematological profiles are associated with CD4-defined disease status in HIV infection.”**

2. **“Systemic inflammation and immune activation in HIV infection: IL-6 is associated with CD4-defined immunodeficiency and sIL-2R declines early after antiretroviral therapy initiation.”**

The purpose of this repository is to ensure **transparency, reproducibility, and reuse** of the analyses performed in these studies.

---

## 📌 Study overview

Human immunodeficiency virus (HIV) infection is characterized by chronic immune activation and systemic inflammation, which contribute to progressive immunodeficiency and non-AIDS–related comorbidities.  
This study combines **cross-sectional** and **longitudinal** analyses to investigate soluble immune mediators and classical immunological markers in people living with HIV (PLWH).

The analyses include:
- Immune stratification by CD4⁺ T-cell count
- Exploratory multivariate analyses (PCA and K-means clustering)
- Multivariate logistic regression to identify predictors of advanced immunodeficiency
- Longitudinal paired analyses evaluating early immune modulation after antiretroviral therapy (ART)

All patients included in the longitudinal analysis received a **dolutegravir-based ART regimen (tenofovir + lamivudine + dolutegravir)**.

---

## 📂 Repository structure
```
HIV-immunological-evaluation/
├── scripts/
│ ├── data_cleaning.R
│ ├── cd4_group_analysis.R
│ ├── pca_clustering.R
│ ├── logistic_regression.R
│ └── treatment_longitudinal_analysis.R
│
├── data/
│ ├── hiv_baseline_data.csv
│ ├── hiv_treatment_data.csv
│ └── README.md
│
├── LICENSE
├── LICENSE_DATA
└── README.md
```

## 🧪 Data

The `data/` directory contains anonymized datasets used in the analyses.

- All data were **fully anonymized** prior to sharing.
- No direct or indirect personal identifiers are included.
- Data are provided exclusively for **research and reproducibility purposes**.

See `data/README.md` for details on variables and structure.

---

## 📊 Analyses included

The scripts reproduce all analyses reported in the manuscript, including:

- Descriptive and comparative analyses across CD4⁺ T-cell strata
- Non-parametric statistical testing
- Principal Component Analysis (PCA)
- K-means and hierarchical clustering
- Multivariate logistic regression with univariate screening (p ≤ 0.20)
- Longitudinal paired analyses (Wilcoxon signed-rank test)
- Generation of publication-ready figures and tables

All scripts were written in **R (version ≥ 4.3)**.

---

## 🔁 Reproducibility

To reproduce the analyses:

1. Clone the repository:
   ```bash
   git clone https://github.com/Bruno-4lmeida/HIV-immunological-evaluation.git
   ```
---
## 📜 License

- Code (scripts directory): MIT License  
- Data (data directory): Creative Commons Attribution 4.0 International (CC BY 4.0)

See `LICENSE` and `LICENSE_DATA` for full license texts.


## 📖 Citation

If you use this repository, please cite the associated publication and this repository.


## 👤 Author

Bruno Almeida Silva
Biomedical Scientist | Immunology | Bioinformatics
Federal University of Pernambuco (UFPE)
IRCCS Burlo Garofolo
