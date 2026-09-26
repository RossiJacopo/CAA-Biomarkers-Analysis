# Cerebral Amyloid Angiopathy (CAA) - Predictive Modeling Pipeline

This repository contains advanced statistical and Machine Learning pipelines developed in R to analyze clinical and radiological biomarkers for Cerebral Amyloid Angiopathy (CAA). 

The project evaluates both Cerebrospinal Fluid (CSF) and Plasma biomarkers, utilizing robust predictive modeling, stability selection, and ensemble learning techniques to handle limited sample sizes and prevent overfitting.

## 📂 Repository Structure

The workflow is divided into three main scripts:

*   **`CAA_Statistical_and_Predictive_Pipeline.r`** 
    Contains the full exploratory data analysis (EDA), multi-variate statistical testing (MANOVA with Bonferroni confidence intervals), and the baseline logistic regression models targeting radiological definitions[cite: 7]. Includes strict implementation of Yeo-Johnson transformations and cross-validation[cite: 7].
*   **`CAA_Predictive_Model_Lumipulse.r`** 
    Focuses on Lumipulse biomarkers. It implements a complete predictive pipeline featuring robust data augmentation (Synthetic Minority Oversampling via convex combinations and KNN), variable selection (Lasso, Exhaustive Search), and a Generalized Linear Mixed Model (GLMM) handling random effects[cite: 5].
*   **`CAA_Stability_Selection_and_Bagging.r`** 
    Dedicated to Plasma biomarkers and ensemble methods. It features a Stability Selection process running Lasso regression across 100 bootstrap samples to find the most robust predictors[cite: 6]. The final predictions are made using an Ensemble Bagging approach (100 bagged GLM models) evaluated through Repeated Cross-Validation and LOOCV[cite: 6].

## 🛠️ Methodological Highlights
*   **Data Leakage Prevention:** Synthetic data augmentation and feature scaling are strictly confined within the training sets during data partitioning and cross-validation folds.
*   **Robust Feature Selection:** Instead of relying on a single penalized regression, variables are selected based on their survival frequency across bootstrapped Lasso iterations (>80% selection rate).
*   **Dynamic Thresholding:** Classification thresholds are dynamically calculated using Youden's Index to maximize Sensitivity and Specificity.

*Note: The original raw clinical datasets (`dataset_1.xlsx`, `data2.RData`) are not included in this repository to comply with GDPR and medical data privacy policies.*
