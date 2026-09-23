# =====================================================================
# SCRIPT PROGETTO.R 
# =====================================================================

# ---- Load packages ----
library(mvtnorm); library(car); library(readxl); library(ellipse)
library(faraway); library(leaps); library(MASS); library(GGally)
library(rgl); library(dplyr); library(RColorBrewer); library(MVN)
library(naniar); library(caret); library(pROC); library(glmnet)
library(ggplot2); library(insight); library(lattice); library(lme4)
library(arm); library(ResourceSelection); library(heplots)

# ---- Helper Functions ----

# Plot Histograms and QQ-plots
plot_hist_qq <- function(df) {
  n_cols <- ncol(df)
  par(mfrow = c(n_cols, 2))
  for (i in 1:n_cols) {
    x <- df[, i]
    name <- colnames(df)[i]
    hist(x, prob = TRUE, col = 'grey85', main = paste('Histogram', name))
    lines(seq(min(x), max(x), length.out = 100), dnorm(seq(min(x), max(x), length.out = 100), mean(x), sd(x)), col = 'blue', lty = 2)
    qqnorm(x, main = paste('QQplot', name))
    qqline(x)
  }
  par(mfrow = c(1, 1))
}

# Plot Covariance Heatmaps
plot_cov_matrices <- function(data_num, group_var) {
  cov_list <- lapply(split(data_num, group_var), cov)
  all_matrices <- do.call(rbind, cov_list)
  breaks <- quantile(all_matrices, (0:100)/100, na.rm = TRUE)
  par(mfrow = c(1, length(cov_list)))
  for(i in seq_along(cov_list)) {
    image(cov_list[[i]], col = heat.colors(100), main = paste('Cov.', names(cov_list)[i]), asp = 1, axes = FALSE, breaks = breaks)
  }
  par(mfrow = c(1, 1))
}

# Compute and Plot Bonferroni CIs
compute_and_plot_bonferroni <- function(fit, data_num, group_var, alpha = 0.05) {
  groups <- split(data_num, group_var)
  ns <- sapply(groups, nrow)
  n <- sum(ns); g <- length(groups); p <- ncol(data_num)
  k <- p * g * (g - 1) / 2
  qT <- qt(1 - alpha / (2 * k), n - g)
  
  W <- summary.manova(fit)$SS$Residuals
  means <- sapply(groups, colMeans) 
  
  group_names <- names(groups)
  pairs <- combn(group_names, 2, simplify = FALSE)
  pair_names <- sapply(pairs, paste, collapse = "-")
  
  par(mfrow = c(1, p), mar = c(8, 4, 4, 2))
  for (j in 1:p) {
    lowers <- numeric(length(pairs)); uppers <- numeric(length(pairs)); diffs <- numeric(length(pairs))
    for (i in seq_along(pairs)) {
      g1 <- pairs[[i]][1]; g2 <- pairs[[i]][2]
      diff_val <- means[j, g1] - means[j, g2]
      se <- sqrt(W[j, j] / (n - g) * (1/ns[g1] + 1/ns[g2]))
      diffs[i] <- diff_val; lowers[i] <- diff_val - qT * se; uppers[i] <- diff_val + qT * se
    }
    ylim <- c(min(lowers) - 0.1 * diff(range(lowers, uppers)), max(uppers) + 0.1 * diff(range(lowers, uppers)))
    plot(1:length(pairs), ylim = ylim, xlim = c(0.5, length(pairs) + 0.5), pch = '', xlab = '', ylab = paste('CI', colnames(data_num)[j]), main = paste('CI', colnames(data_num)[j]), xaxt = "n")
    axis(1, at = 1:length(pairs), labels = FALSE)
    text(x = 1:length(pairs), y = par("usr")[3] - 0.05 * diff(par("usr")[3:4]), labels = pair_names, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
    for (i in seq_along(pairs)) {
      lines(c(i, i), c(lowers[i], uppers[i])); points(i, diffs[i], pch = 16)
      points(i, lowers[i], col = "blue", pch = 16); points(i, uppers[i], col = "red", pch = 16)
    }
    abline(h = 0, lty = 2)
  }
  par(mfrow = c(1, 1))
}

# Synthetic data augmentation
augment_data <- function(data, target_col = "target", numeric_cols = NULL, categorical_cols = NULL, n_new = 20, k_range = c(2, 5), noise_sd = 0.05, k_nn = 3) {
  class0 <- data[data[[target_col]] == 0, ]
  n_class0 <- nrow(class0)
  if (n_class0 < max(k_range)) stop("Not enough observations in class 0.")
  if (is.null(numeric_cols)) stop("Please specify numeric columns.")
  
  augmented <- data.frame()
  for (i in 1:n_new) {
    k <- sample(k_range[1]:k_range[2], 1)
    idx <- sample(1:n_class0, k, replace = FALSE)
    selected_numeric <- class0[idx, numeric_cols, drop = FALSE]
    weights <- runif(k); weights <- weights / sum(weights)
    combined_numeric <- colSums(selected_numeric * weights)
    noisy_combined <- combined_numeric + rnorm(length(combined_numeric), mean = 0, sd = noise_sd)
    
    new_point <- as.list(noisy_combined)
    if (!is.null(categorical_cols) && length(categorical_cols) > 0) {
      distances <- apply(class0[, numeric_cols, drop = FALSE], 1, function(row) sqrt(sum((row - noisy_combined)^2)))
      nearest_indices <- order(distances)[1:k_nn]
      nearest_cat <- class0[nearest_indices, categorical_cols, drop = FALSE]
      assigned_cat <- sapply(nearest_cat, function(col) names(which.max(table(col))))
      new_point <- c(new_point, as.list(assigned_cat))
    }
    new_point[[target_col]] <- 0
    new_df <- as.data.frame(new_point, stringsAsFactors = FALSE)
    for (v in numeric_cols) new_df[[v]] <- as.numeric(new_df[[v]])
    if (!is.null(categorical_cols)) for (v in categorical_cols) new_df[[v]] <- type.convert(new_df[[v]], as.is = TRUE)
    augmented <- rbind(augmented, new_df)
  }
  return(augmented[, colnames(data)])
}

# ---- Read dataset ----
data <- read_excel("dataset_1.xlsx")
data_CAA <- data[-(192:237), ]

# ---- Preliminary analysis ----
cols_to_remove <- c(
  "Date of birth", "Year of admission (outpatient clinic or hospitalization) in our Institute",
  "T0 (year)_VAL", "T1 (year)_VAL", "T2 (year)_VAL", "Age at admission in our institute_RECRUIT", 
  "Age at recruitement in this study_VAL", "antiplatelet FU (type)", "MTA_RM_REC", "PTA_RM_REC", 
  "CST3_GEN", "ITM2B_GEN", "TTR_GEN", "z-score_PET", "Event (3)_ev_clinica_ffe", "Year 3_ev_clinica_ffe",
  "Event (4)_ev_clinica_ffe", "Year 4_ev_clinica_ffe", "Event (5)_ev_clinica_ffe", "Year 5_ev_clinica_ffe",
  "Event (6)_ev_clinica_ffe", "Year 6_ev_clinica_ffe", "Duration retrospective and prospective FU (months)",
  "NfL (pg/ml)CFST1", "GFAP(pg/ml)CFST1", "T0_sdt_DISABILITA", "T1_sdt_DISABILITA", "T2_sdt_DISABILITA",
  "data_mocaT0_SCOREDEC", "MoCA t2_SCOREDEC", "P.G._SCOREDEC...34", "P.G._SCOREDEC...26", 
  "P.G._SCOREDEC...30", "P.C._SCOREDEC...35", "P.E._SCOREDEC...36", "P.E._SCOREDEC...32", "P.E._SCOREDEC...28",
  "Last MRI available (year)_RM_REC", "1st MRI available (year)_RM_REC", "Year 1_ev_clinica_ffe", 
  "Year 2_ev_clinica_ffe", "Data esecuzio PL_CFSlumipulse", "MRI available (number)_RM_REC",
  "CSF Dinamica / AD like_RM_REC", "CAA type (i-CAA / CAA-ri / CAA / h-CAA)_CD", "Possible 2.0_CD",
  "Probable 2.0_CD", "sICH, TF, cognitive impairment/dementia_CD", "MRI-proven ICH, CMB, cSSS, cSAH_CD",
  "MRI-proven lobar lesion+white matter lesions (WMH or PVS)_CD", "APP_GEN", "other_GEN", "PET-imaging",
  "tracer_PET", "NfL (pg/ml)CFS", "GFAP(pg/ml)CFS", "Ethnic group", "MoCA t0_SCOREDEC", "MoCA t1_SCOREDEC",
  "any SAH_RM_REC", "ApoE4 carrier_GEN", "ApoE2 carrier_GEN"
)
data2 <- data_CAA %>% select(-one_of(cols_to_remove))

data2$`antiplatelet (type)_TERAPIA_REC` <- as.numeric(gsub(" ", "", gsub("%", "", gsub(",", ".", data2$`antiplatelet (type)_TERAPIA_REC`))))
data2$`anticoagulants (type)_TERAPIA_REC` <- as.numeric(data2$`anticoagulants (type)_TERAPIA_REC`)
data2$Tesla_RM_REC <- as.numeric(data2$Tesla_RM_REC) data2$`Aβ42/Aβ40_CFSlumipulse` <- as.numeric(data2$`Aβ42/Aβ40_CFSlumipulse`)
data2$`Aβ40(pg/ml)PLASMA` <- as.numeric(data2$`Aβ40(pg/ml)PLASMA`)
data2$`NfL (pg/ml)PLASMA` <- as.numeric(data2$`NfL (pg/ml)PLASMA`)
data2$Genetics_GEN <- as.numeric(data2$Genetics_GEN)

factor_cols <- c(
  "ICH predominance side_RM_REC", "Index event type", "MRS before index event_DISABILITA",
  "MRS at T0 (1a valutazio al Besta)_DISABILITA", "MRS at T1_DISABILITA", "MRS at T2_DISABILITA",
  "AF (0=no; 1=yes; 2=AF after diagnosis)_RISK", "Smoking (0=no;1=yes;2=former)_RISK",
  "Alcohol(0=no;1=yes;2=former_RISK", "lobar CMB (0-3)_RM_REC", "deep CMB (0-2)_RM_REC",
  "cSS_RM_REC", "WMHs_RM_REC", "Confluent WMHs_RM_REC", "Fazekas periventricular WM_RM_REC",
  "Fazekas deep WM_RM_REC", "Scheltens' scale_RM_REC", "Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC", 
  "SWI / FFE-GRE_RM_REC", "CAA type (familial/sporadic)_CD", "ApoEallele1_GEN", "ApoEallele2_GEN", 
  "Event (1)_ev_clinica_ffe", "Event (2)_ev_clinica_ffe", "ICH (1/0)_RECRUIT", "Cognitive impairment (1/0)_RECRUIT",
  "Dementia (1/0)_RECRUIT", "TF (0=no;1=yes)_RECRUIT", "Sex (0=male;1=female)", "Hypertension_RISK", 
  "Dyslipidemia_RISK", "Diabetes_RISK", "Auricola_RISK", "Previous stroke_RISK", "Previous ICH(0=no;1=yes)_RISK", 
  "Autoimmu disorders(0=no;1=yes)_RISK", "previous neurosurgery", "antiplatelet (type)_TERAPIA_REC", 
  "anticoagulants (type)_TERAPIA_REC", "statins(0=no;1=yes)_TERAPIA_REC", "antiseizure medications(0=no;1=yes)_TERAPIA_REC", 
  "lobar ICH_RM_REC", "Cerebellar ICH_RM_REC", "Cerebellar SS_RM_REC", "Cerebellar lacu / stroke_RM_REC", 
  "Posterior predominance_RM_REC", "Multispot WM hyperintensity pattern (>10)_RM_REC", 
  "Centrum Semiovale-PVS (>20/emisfero)_RM_REC", "DEEP-PVS (basal ganglia, pons)_RM_REC", 
  "convexity SAH_RM_REC", "DWI+ lesions_RM_REC", "deep lacunae (1/0)_RM_REC", "Genetics_GEN", 
  "result_PET", "Cognitive detrimental during FU"
)

for (col in factor_cols) {
  if (col %in% colnames(data2)) data2[[col]] <- as.factor(data2[[col]])
}

numeric_columns <- data2[sapply(data2, is.numeric)]
numeric_columns <- numeric_columns[, -8]

# ---- PCA Analysis ----
numeric_bio <- na.omit(numeric_columns[, c(8,9,11,12)])
PCA_numeric_bio_std <- princomp(scale(numeric_bio))

par(mfrow = c(1, 2))
plot(PCA_numeric_bio_std, las = 2, main = 'Principal Components', ylim = c(0, 7), col = 'orange', lwd = 2)
plot(cumsum(PCA_numeric_bio_std$sde^2) / sum(PCA_numeric_bio_std$sde^2), type = 'b', axes = FALSE, 
     xlab = 'Number of components', ylab = 'Contribution to the total variance', ylim = c(0, 1))
abline(h = 1, col = 'forestgreen')
box(); axis(2); axis(1)
par(mfrow = c(1, 1))

# ---- SECTION 1: MANOVA BIO MARKERS & DISABILITY ----
manova_bio <- na.omit(cbind(data2[, 10], numeric_columns[, c(8,9,11,12)]))
manova_bio$`MRS at T0 (1a valutazio al Besta)_DISABILITA` <- as.factor(manova_bio$`MRS at T0 (1a valutazio al Besta)_DISABILITA`)

b <- manova_bio[, 2:5]
group <- manova_bio[, 1]
M <- colMeans(b); S <- cov(b)
d2 <- matrix(mahalanobis(b, M, S))
i_valid <- which(d2 <= 7)

b_wo_outliers <- b[i_valid, ]
group_wo <- group[i_valid]

lambda <- powerTransform(b_wo_outliers)
b_transformed <- b_wo_outliers
b_transformed[, 1] <- bcPower(b_wo_outliers[, 1], lambda$lambda[1])
b_transformed[, 2] <- bcPower(b_wo_outliers[, 2], lambda$lambda[2])
b_transformed[, 3] <- yjPower(b_wo_outliers[, 3], powerTransform(b_wo_outliers[, 3], family = "yjPower")$lambda)
b_transformed[, 4] <- yjPower(b_wo_outliers[, 4], powerTransform(b_wo_outliers[, 4], family = "yjPower")$lambda)

plot_hist_qq(b_transformed)
plot_cov_matrices(b_transformed, group_wo)

fit_bio <- manova(as.matrix(b_transformed) ~ group_wo)
summary.aov(fit_bio)
compute_and_plot_bonferroni(fit_bio, b_transformed, group_wo)

# ---- SECTION 2: MANOVA RADIOLOGICAL DEFINITION (PLASMA) ----
manova_CAA <- na.omit(data.frame(
  group = as.factor(data_CAA$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC`),
  Ab40 = data_CAA$`Aβ40(pg/ml)PLASMA`, 
  Ab42 = data_CAA$`Aβ42 (pg/ml)PLASMA`, 
  Tau = data_CAA$`Tau (pg/ml)PLASMA`
))

b <- manova_CAA[, 2:4]
group <- manova_CAA[, 1]
M <- colMeans(b); S <- cov(b)
d2 <- matrix(mahalanobis(b, M, S))
i_valid <- which(d2 <= 9)

b_wo_outliers <- b[i_valid, ]
group_wo <- group[i_valid]

lambda <- powerTransform(b_wo_outliers)
b_transformed <- b_wo_outliers
for (i in 1:3) b_transformed[, i] <- bcPower(b_wo_outliers[, i], lambda$lambda[i])

plot_hist_qq(b_transformed)
plot_cov_matrices(b_transformed, group_wo)

fit_caa <- manova(as.matrix(b_transformed) ~ group_wo)
summary.aov(fit_caa)
compute_and_plot_bonferroni(fit_caa, b_transformed, group_wo)

# ---- SECTION 3: MANOVA RADIOLOGICAL DEFINITION (CFS) ----
manova_CFS <- na.omit(data.frame(
  group = as.factor(data2$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC`),
  Ab40 = data2$`Aβ40(pg/ml)CFS`,
  Ab42 = data2$`Aβ42 (pg/ml)CFS`, 
  Tau = data2$`Tau (pg/ml)CFS`
))

b <- manova_CFS[, 2:4]
group <- manova_CFS[, 1]
M <- colMeans(b); S <- cov(b)
d2 <- matrix(mahalanobis(b, M, S))
i_valid <- which(d2 <= 8)

b_wo_outliers <- b[i_valid, ]
group_wo <- group[i_valid]

lambda <- powerTransform(b_wo_outliers)
b_transformed <- b_wo_outliers
for (i in 1:3) b_transformed[, i] <- bcPower(b_wo_outliers[, i], lambda$lambda[i])

plot_hist_qq(b_transformed)
plot_cov_matrices(b_transformed, group_wo)

fit_cfs <- manova(as.matrix(b_transformed) ~ group_wo)
summary.aov(fit_cfs)
compute_and_plot_bonferroni(fit_cfs, b_transformed, group_wo)

# ---- SECTION 4: BI-MANOVA RADIOLOGICAL DEFINITION & SEX ----
bi_manova_lumipulse <- na.omit(data.frame(
  rad = factor(ifelse(data2$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC` %in% c("CAA", "CAA-DPA", "iCAA"), 1, 0), labels = c("Non CAA", "CAA")),
  sex = as.factor(data2$`Sex (0=male;1=female)`),
  age = data2$`Age at onset of CAA/DPA symptoms_RECRUIT`,
  scoredec = data2$P.C._SCOREDEC...27,   Ab40 = data2$`Abeta40 - Valore_CFSlumipulse`, 
  Ab42 = data2$`Abeta42 - Valore_CFSlumipulse`, 
  Tau = data2$`total-tau - Valore_CFSlumipulse`
))

bi_manova_num <- bi_manova_lumipulse[, 3:7]
Mx <- colMeans(bi_manova_num); Sx <- cov(bi_manova_num)
d2x <- matrix(mahalanobis(bi_manova_num, Mx, Sx))
i_validx <- which(d2x <= 13)

bi_manova_lumipulse <- bi_manova_lumipulse[i_validx, ]
bi_manova_num <- bi_manova_lumipulse[, 3:7]

lambdax <- powerTransform(bi_manova_num)
for(i in 1:5) bi_manova_lumipulse[, i+2] <- bcPower(bi_manova_lumipulse[, i+2], lambdax$lambda[i])

bi_manova_num <- scale(bi_manova_lumipulse[, 3:7])
bi_manova_lumipulse[, 3:7] <- bi_manova_num

pat_sex <- factor(paste(bi_manova_lumipulse$rad, bi_manova_lumipulse$sex))
plot_cov_matrices(bi_manova_num, pat_sex)

fit_bi <- manova(as.matrix(bi_manova_num) ~ rad + sex + rad:sex, data = bi_manova_lumipulse)
summary.manova(fit_bi)
fit_bi_add <- manova(as.matrix(bi_manova_num) ~ rad + sex, data = bi_manova_lumipulse)
summary.aov(fit_bi_add)
compute_and_plot_bonferroni(fit_bi_add, bi_manova_num, pat_sex)

# ---- SECTION 5: BASIC LOGISTIC REGRESSION EXPLORATION ----
data_regression <- data.frame(
  target = data2$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC`,
  AB40_lumipulse = data2$`Abeta40 - Valore_CFSlumipulse`,
  AB42_lumipulse = data2$`Abeta42 - Valore_CFSlumipulse`,
  Tau_lumipulse = data2$`p-tau181 - Valore_CFSlumipulse`
)

data_regression$target <- ifelse(data_regression$target %in% c("CAA", "CAA-DPA", "iCAA"), 1, 0)
data_regression$target <- factor(data_regression$target, levels = c(0, 1))
data_regression <- na.omit(data_regression)

data_regression$AB40_lumipulse <- scale(data_regression$AB40_lumipulse)
data_regression$AB42_lumipulse <- scale(data_regression$AB42_lumipulse)
data_regression$Tau_lumipulse <- scale(data_regression$Tau_lumipulse)

numeric_data <- data_regression[, 2:4]
M <- colMeans(numeric_data); S <- cov(numeric_data)
d2 <- matrix(mahalanobis(numeric_data, M, S))
valid_indexes <- which(d2 <= 13)
numeric_data <- numeric_data[valid_indexes, ]
data_regression <- data_regression[valid_indexes, ]

M <- colMeans(numeric_data); S <- cov(numeric_data)
d2 <- matrix(mahalanobis(numeric_data, M, S))
valid_indexes <- which(d2 <= 20)
data_regression <- data_regression[valid_indexes, ]

model1.0 <- glm(target ~ AB40_lumipulse + AB42_lumipulse + Tau_lumipulse, data = data_regression, family = 'binomial')
model1.1 <- glm(target ~ AB42_lumipulse + Tau_lumipulse + Tau_lumipulse:AB42_lumipulse, data = data_regression, family = 'binomial')

classe_0 <- data_regression[data_regression$target == 0, ]            classe_1 <- data_regression[data_regression$target == 1, ] 

set.seed(123)
classe_1_sample_40 <- classe_1[sample(nrow(classe_1), size = 40), ]
balanced_data_reg_40 <- rbind(classe_0, classe_1_sample_40)
model2.0 <- glm(target ~ AB42_lumipulse + Tau_lumipulse + Tau_lumipulse:AB42_lumipulse, data = balanced_data_reg_40, family = 'binomial')

set.seed(123)
classe_1_sample_20 <- classe_1[sample(nrow(classe_1), size = 20), ]
balanced_data_reg_20 <- rbind(classe_0, classe_1_sample_20)
model2.1 <- glm(target ~ AB42_lumipulse + Tau_lumipulse + Tau_lumipulse:AB42_lumipulse, data = balanced_data_reg_20, family = 'binomial')

model3.0 <- glm(target ~ AB40_lumipulse + AB42_lumipulse + Tau_lumipulse + AB40_lumipulse:AB42_lumipulse + AB40_lumipulse:Tau_lumipulse + Tau_lumipulse:AB42_lumipulse, data = balanced_data_reg_20, family = 'binomial')
model3.1 <- glm(target ~ AB40_lumipulse + AB42_lumipulse + Tau_lumipulse + AB40_lumipulse:AB42_lumipulse + AB40_lumipulse:Tau_lumipulse + Tau_lumipulse:AB42_lumipulse, data = balanced_data_reg_40, family = 'binomial')

# ---- SECTION 6: EXTENDED MODEL ("MODELLONE") ----
data_modellone <- na.omit(data.frame(
  target = factor(ifelse(data2$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC` %in% c("CAA", "CAA-DPA", "iCAA"), 1, 0), levels = c(0, 1)),
  AB40_lumipulse = data2$`Abeta40 - Valore_CFSlumipulse`,
  AB42_lumipulse = data2$`Abeta42 - Valore_CFSlumipulse`,
  Tau_lumipulse = data2$`p-tau181 - Valore_CFSlumipulse`,
  Age_onset_symptoms = data2$`Age at onset of CAA/DPA symptoms_RECRUIT`,
  CMBS_deep = data2$`CMBs (deep)_RM_REC`,
  CMBS_lobar = data2$`CMBs (lobar) _RM_REC`
))

numeric_modellone <- data_modellone[, -1]
M <- colMeans(numeric_modellone); S <- cov(numeric_modellone)
d2 <- matrix(mahalanobis(numeric_modellone, M, S))
data_modellone <- data_modellone[which(d2 <= 30), ]

lambdam <- powerTransform(data_modellone[, -1], family = "yjPower")
for (i in 2:ncol(data_modellone)) {
  data_modellone[, i] <- yjPower(data_modellone[, i], lambdam$lambda[i - 1])
}

modellone_1 <- glm(target ~ ., data = scale(data_modellone[, -1]) %>% as.data.frame() %>% mutate(target = data_modellone$target), family = 'binomial')
modellone_2 <- glm(target ~ AB42_lumipulse + Tau_lumipulse + Age_onset_symptoms + CMBS_deep + CMBS_lobar, data = data_modellone, family = 'binomial')
modellone_3 <- glm(target ~ AB42_lumipulse + Tau_lumipulse + CMBS_deep + CMBS_lobar, data = data_modellone, family = 'binomial')
modellone_4 <- glm(target ~ AB42_lumipulse + Tau_lumipulse + CMBS_deep + CMBS_lobar + AB42_lumipulse:Tau_lumipulse + AB42_lumipulse:CMBS_deep, data = data_modellone, family = 'binomial')

set.seed(123)
train_idx <- createDataPartition(data_modellone$target, p = 0.7, list = FALSE)
train_data <- data_modellone[train_idx, ]
test_data <- data_modellone[-train_idx, ]

numeric_vars_mod <- names(train_data)[-1]
set.seed(123)
zeri_augmented_modellone <- augment_data(train_data, target_col = "target", numeric_cols = numeric_vars_mod, categorical_cols = NULL, n_new = 26)
train_data_augmented <- rbind(train_data, zeri_augmented_modellone)
train_data_augmented <- train_data_augmented[sample(nrow(train_data_augmented)), ]

preproc <- preProcess(train_data_augmented[, -1], method = c("center", "scale"))
train_data_augmented[, -1] <- predict(preproc, train_data_augmented[, -1])
test_data[, -1] <- predict(preproc, test_data[, -1])

modellone_4_trained <- glm(target ~ AB42_lumipulse + CMBS_deep + CMBS_lobar +
                             AB42_lumipulse:CMBS_deep + AB42_lumipulse:CMBS_lobar +
                             CMBS_deep:CMBS_lobar, data = train_data_augmented, family = "binomial")

prob_test <- predict(modellone_4_trained, newdata = test_data, type = "response")
roc_obj <- roc(test_data$target, prob_test)
best_thresh <- coords(roc_obj, "best", ret="threshold")$threshold[1]
pred_class_test <- factor(ifelse(prob_test >= best_thresh, 1, 0), levels = c(0, 1))

confusionMatrix(pred_class_test, factor(test_data$target, levels = c(0, 1)), positive = "1")
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve - Modellone Test Set")

data_cv <- rbind(train_data_augmented, test_data)
data_cv$target <- factor(ifelse(data_cv$target == 1, "X1", "X0"))
ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")

modellone_cv <- train(
  target ~ AB42_lumipulse + CMBS_deep + CMBS_lobar + AB42_lumipulse:CMBS_deep + AB42_lumipulse:CMBS_lobar + CMBS_deep:CMBS_lobar,
  data = data_cv, method = "glm", family = "binomial", 
  preProcess = c("center", "scale"), trControl = ctrl, metric = "ROC"
)

# ---- SECTION 7: PLASMA MODEL WITH CATEGORICALS ("MODELLISSIMO") ----
data_modellissimo <- data.frame(
  target = data2$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC`,
  AB40_plasma = data2$`Aβ40(pg/ml)PLASMA`,
  AB42_plasma = data2$`Aβ42 (pg/ml)PLASMA`,
  Tau_plasma = data2$`Tau (pg/ml)PLASMA`,
  Age_onset_symptoms = data2$`Age at onset of CAA/DPA symptoms_RECRUIT`,
  CMBS_deep = data2$`CMBs (deep)_RM_REC`,
  CMBS_lobar = data2$`CMBs (lobar) _RM_REC`,
  Sex = data2$`Sex (0=male;1=female)`,
  Dyslipidemia_RISK = data2$Dyslipidemia_RISK,   Neurosurgery = data2$`previous neurosurgery`,
  Hypertension = data2$Hypertension_RISK,   ICH = data2$`ICH (1/0)_RECRUIT`,
  Stroke = data2$`Previous stroke_RISK`,
  Alchool = data2$`Alcohol(0=no;1=yes;2=former_RISK`,
  AF = data2$`AF (0=no; 1=yes; 2=AF after diagnosis)_RISK`,
  Smoking = data2$`Smoking (0=no;1=yes;2=former)_RISK`,
  Antiplatelet = data2$`antiplatelet (type)_TERAPIA_REC`,
  Anticoagulants = data2$`anticoagulants (type)_TERAPIA_REC`
)

data_modellissimo$target <- ifelse(data_modellissimo$target %in% c("CAA", "CAA-DPA", "iCAA"), 1, 0)
data_modellissimo$target <- factor(data_modellissimo$target, levels = c(0, 1))
data_modellissimo <- na.omit(data_modellissimo)

numeric_modellissimo <- data_modellissimo[, 2:7]
M <- colMeans(numeric_modellissimo); S <- cov(numeric_modellissimo)
d2 <- matrix(mahalanobis(numeric_modellissimo, M, S))
data_modellissimo <- data_modellissimo[which(d2 <= 20), ]

regfit.fwd <- regsubsets(target ~ ., data = data_modellissimo, nvmax = 11, method = "forward")
modellissimo_1 <- glm(target ~ CMBS_deep + AB42_plasma + Alchool + Anticoagulants + CMBS_lobar + Dyslipidemia_RISK + AB40_plasma + AF + Stroke + ICH, data = data_modellissimo, family = 'binomial')

forced_vars <- c("AB40_plasma", "AB42_plasma", "CMBS_deep", "CMBS_lobar")
predictor_names <- names(data_modellissimo)[names(data_modellissimo) != "target"]
forced_indices <- which(predictor_names %in% forced_vars)
regfit.full <- regsubsets(target ~ ., data = data_modellissimo, nvmax = 11, force.in = forced_indices)

vars_model <- c("target", "AB42_plasma", "CMBS_deep", "CMBS_lobar", "ICH")
data_mod_sub <- data_modellissimo[, vars_model]

set.seed(123)
train_idx_mod <- createDataPartition(data_mod_sub$target, p = 0.7, list = FALSE)
train_data_mod <- data_mod_sub[train_idx_mod, ]
test_data_mod <- data_mod_sub[-train_idx_mod, ]

numeric_vars_mod <- c("AB42_plasma", "CMBS_deep", "CMBS_lobar")
categorical_vars_mod <- c("ICH")

set.seed(123)
zeri_aug_mod <- augment_data(train_data_mod, target_col = "target", numeric_cols = numeric_vars_mod, categorical_cols = categorical_vars_mod, n_new = 20)
train_data_aug_mod <- rbind(train_data_mod, zeri_aug_mod)
train_data_aug_mod <- train_data_aug_mod[sample(nrow(train_data_aug_mod)), ]

lambdam_yj <- powerTransform(train_data_aug_mod[, c("CMBS_deep", "CMBS_lobar")], family = "yjPower")

train_data_aug_mod$CMBS_deep_yj <- yjPower(train_data_aug_mod$CMBS_deep, lambdam_yj$lambda[1])
train_data_aug_mod$CMBS_lobar_yj <- yjPower(train_data_aug_mod$CMBS_lobar, lambdam_yj$lambda[2])

test_data_mod$CMBS_deep_yj <- yjPower(test_data_mod$CMBS_deep, lambdam_yj$lambda[1])
test_data_mod$CMBS_lobar_yj <- yjPower(test_data_mod$CMBS_lobar, lambdam_yj$lambda[2])

preproc_mod <- preProcess(train_data_aug_mod[, c("AB42_plasma", "CMBS_deep_yj", "CMBS_lobar_yj")], method = c("center", "scale"))
train_data_aug_mod[, c("AB42_plasma", "CMBS_deep_yj", "CMBS_lobar_yj")] <- predict(preproc_mod, train_data_aug_mod[, c("AB42_plasma", "CMBS_deep_yj", "CMBS_lobar_yj")])
test_data_mod[, c("AB42_plasma", "CMBS_deep_yj", "CMBS_lobar_yj")] <- predict(preproc_mod, test_data_mod[, c("AB42_plasma", "CMBS_deep_yj", "CMBS_lobar_yj")])

modellissimo_YJ <- glm(target ~ AB42_plasma + CMBS_deep_yj + ICH + 
                         AB42_plasma:CMBS_deep_yj + AB42_plasma:CMBS_lobar_yj,
                       data = train_data_aug_mod, family = 'binomial')

prob_test_mod <- predict(modellissimo_YJ, newdata = test_data_mod, type = "response")
roc_obj_mod <- roc(test_data_mod$target, prob_test_mod)
best_thresh_mod <- coords(roc_obj_mod, "best", ret="threshold")$threshold[1]
pred_class_mod <- factor(ifelse(prob_test_mod >= best_thresh_mod, 1, 0), levels = c(0, 1))

conf_obj_mod <- confusionMatrix(pred_class_mod, factor(test_data_mod$target, levels = c(0, 1)), positive = '1')

plot(roc_obj_mod, col = "blue", main = "ROC Curve - Modellissimo Test Set")

binnedplot(fitted(modellissimo_YJ), rstandard(modellissimo_YJ), main = "Binned Residual Plot", xlab = "Predicted Values", ylab = "Standardized Residuals")

train_data_aug_mod$target_factor <- factor(ifelse(train_data_aug_mod$target == 1, "Yes", "No"), levels = c("Yes", "No"))
formula_cv <- target_factor ~ AB42_plasma + CMBS_deep_yj + ICH + AB42_plasma:CMBS_deep_yj + AB42_plasma:CMBS_lobar_yj

train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
set.seed(123)
cv_model_mod <- train(formula_cv, data = train_data_aug_mod, method = "glm", family = "binomial", trControl = train_control, metric = "ROC")

cm_cv_mod <- confusionMatrix(data = factor(cv_model_mod$pred$pred, levels = c("No", "Yes")),
                             reference = factor(cv_model_mod$pred$obs, levels = c("No", "Yes")), positive = "Yes")

conf_df_cv_mod <- as.data.frame(as.table(cm_cv_mod$table))
colnames(conf_df_cv_mod) <- c("Prediction", "Reference", "Freq")

ggplot(conf_df_cv_mod, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "black") + geom_text(aes(label = Freq), size = 6, color = "white") +
  scale_fill_gradient(low = "pink", high = "red") +
  labs(title = paste0("Confusion Matrix on 5-fold CV\nSensitivity: ", round(cm_cv_mod$byClass["Sensitivity"], 4), " \vert{} Specificity: ", round(cm_cv_mod$byClass["Specificity"], 4)), x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 14)