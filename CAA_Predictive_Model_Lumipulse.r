# ---- Load packages ----
library(mvtnorm)
library(car)
library(readxl)
library(ellipse)
library(faraway)
library(leaps)
library(MASS)
library(GGally)
library(rgl)
library(dplyr)
library(RColorBrewer)
library(MVN)
library(naniar)
library(caret)
library(pROC)
library(glmnet)
library(ggplot2)
library(insight)
library(lattice)
library(lme4)

# ---- Helper Functions ----
# Unified function for synthetic data augmentation
augment_data <- function(data, target_col = "target", numeric_cols = NULL, categorical_cols = NULL, n_new = 20, k_range = c(2, 5), noise_sd = 0.05, k_nn = 3) {
  class0 <- data[data[[target_col]] == 0, ]
  n_class0 <- nrow(class0)
  
  if (n_class0 < max(k_range)) stop("Not enough observations in class 0.")
  if (is.null(numeric_cols)) stop("Please specify numeric columns.")
  
  augmented <- data.frame()
  
  for (i in 1:n_new) {
    # Convex combination for numeric features
    k <- sample(k_range[1]:k_range[2], 1)
    idx <- sample(1:n_class0, k, replace = FALSE)
    selected_numeric <- class0[idx, numeric_cols, drop = FALSE]
    
    weights <- runif(k)
    weights <- weights / sum(weights)
    combined_numeric <- colSums(selected_numeric * weights)
    noisy_combined <- combined_numeric + rnorm(length(combined_numeric), mean = 0, sd = noise_sd)
    
    new_point <- as.list(noisy_combined)
    
    # KNN logic for categorical features
    if (!is.null(categorical_cols) && length(categorical_cols) > 0) {
      distances <- apply(class0[, numeric_cols, drop = FALSE], 1, function(row) {
        sqrt(sum((row - noisy_combined)^2))
      })
      nearest_indices <- order(distances)[1:k_nn]
      nearest_cat <- class0[nearest_indices, categorical_cols, drop = FALSE]
      
      assigned_cat <- sapply(nearest_cat, function(col) names(which.max(table(col))))
      new_point <- c(new_point, as.list(assigned_cat))
    }
    
    new_point[[target_col]] <- 0
    new_df <- as.data.frame(new_point, stringsAsFactors = FALSE)
    
    # Cast variables back to their original types
    for (v in numeric_cols) new_df[[v]] <- as.numeric(new_df[[v]])
    if (!is.null(categorical_cols)) {
      for (v in categorical_cols) new_df[[v]] <- type.convert(new_df[[v]], as.is = TRUE)
    }
    
    augmented <- rbind(augmented, new_df)
  }
  
  # Return preserving original column order
  return(augmented[, colnames(data)])
}

# ---- Read dataset ----
data <- read_excel("dataset_1.xlsx")

# ---- Preliminary analysis ----
data_CAA <- data[-(192:237), ]

cols_to_remove <- c(
  "Date of birth", "Year of admission (outpatient clinic or hospitalization) in our Institute",
  "T0 (year)_VAL", "T1 (year)_VAL", "T2 (year)_VAL", "Age at admission in our institute_RECRUIT",
  "Age at recruitement in this study_VAL", "antiplatelet FU (type)", "MTA_RM_REC", "PTA_RM_REC",
  "CST3_GEN", "ITM2B_GEN", "TTR_GEN", "z-score_PET", "Event (3)_ev_clinica_ffe", "Year 3_ev_clinica_ffe",
  "Event (4)_ev_clinica_ffe", "Year 4_ev_clinica_ffe", "Event (5)_ev_clinica_ffe", "Year 5_ev_clinica_ffe",
  "Event (6)_ev_clinica_ffe", "Year 6_ev_clinica_ffe", "Duration retrospective and prospective FU (months)",
  "NfL (pg/ml)CFST1", "GFAP(pg/ml)CFST1", "T0_sdt_DISABILITA", "T1_sdt_DISABILITA", "T2_sdt_DISABILITA",
  "data_mocaT0_SCOREDEC", "MoCA t2_SCOREDEC", "P.G._SCOREDEC...34", "P.G._SCOREDEC...26", "P.G._SCOREDEC...30",
  "P.C._SCOREDEC...35", "P.E._SCOREDEC...36", "P.E._SCOREDEC...32", "P.E._SCOREDEC...28",
  "Last MRI available (year)_RM_REC", "1st MRI available (year)_RM_REC", "Year 1_ev_clinica_ffe",
  "Year 2_ev_clinica_ffe", "Data esecuzio PL_CFSlumipulse", "MRI available (number)_RM_REC",
  "CSF Dinamica / AD like_RM_REC", "CAA type (i-CAA / CAA-ri / CAA / h-CAA)_CD", "Possible 2.0_CD",
  "Probable 2.0_CD", "sICH, TF, cognitive impairment/dementia_CD", "MRI-proven ICH, CMB, cSSS, cSAH_CD",
  "MRI-proven lobar lesion+white matter lesions (WMH or PVS)_CD", "APP_GEN", "other_GEN", "PET-imaging",
  "tracer_PET", "NfL (pg/ml)CFS", "GFAP(pg/ml)CFS", "Ethnic group", "MoCA t0_SCOREDEC", "MoCA t1_SCOREDEC",
  "any SAH_RM_REC", "ApoE4 carrier_GEN", "ApoE2 carrier_GEN"
)

data2 <- data_CAA %>% select(-one_of(cols_to_remove))

# ---- Data Cleaning & Formatting ----
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
  "Fazekas deep WM_RM_REC", "Scheltens' scale_RM_REC",
  "Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC", "SWI / FFE-GRE_RM_REC",
  "CAA type (familial/sporadic)_CD", "ApoEallele1_GEN", "ApoEallele2_GEN", "Event (1)_ev_clinica_ffe",
  "Event (2)_ev_clinica_ffe", "ICH (1/0)_RECRUIT", "Cognitive impairment (1/0)_RECRUIT",
  "Dementia (1/0)_RECRUIT", "TF (0=no;1=yes)_RECRUIT", "Sex (0=male;1=female)",
  "Hypertension_RISK", "Dyslipidemia_RISK", "Diabetes_RISK", "Auricola_RISK", "Previous stroke_RISK",
  "Previous ICH(0=no;1=yes)_RISK", "Autoimmu disorders(0=no;1=yes)_RISK", "previous neurosurgery",
  "antiplatelet (type)_TERAPIA_REC", "anticoagulants (type)_TERAPIA_REC", "statins(0=no;1=yes)_TERAPIA_REC",
  "antiseizure medications(0=no;1=yes)_TERAPIA_REC", "lobar ICH_RM_REC", "Cerebellar ICH_RM_REC",
  "Cerebellar SS_RM_REC", "Cerebellar lacu / stroke_RM_REC", "Posterior predominance_RM_REC",
  "Multispot WM hyperintensity pattern (>10)_RM_REC", "Centrum Semiovale-PVS (>20/emisfero)_RM_REC",
  "DEEP-PVS (basal ganglia, pons)_RM_REC", "convexity SAH_RM_REC", "DWI+ lesions_RM_REC",
  "deep lacunae (1/0)_RM_REC", "Genetics_GEN", "result_PET", "Cognitive detrimental during FU"
)

for (col in factor_cols) {
  if (col %in% colnames(data2)) data2[[col]] <- as.factor(data2[[col]])
}

# ---- Lumipulse regression model dataset ----
lumipulse_regression <- data.frame(
  target = data2$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC`,
  AB40_lumipulse = data2$`Abeta40 - Valore_CFSlumipulse`,
  AB42_lumipulse = data2$`Abeta42 - Valore_CFSlumipulse`,
  Tau_lumipulse = data2$`p-tau181 - Valore_CFSlumipulse`,
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

lumipulse_regression$target <- ifelse(lumipulse_regression$target %in% c("CAA", "CAA-DPA", "iCAA"), 1, 0)
lumipulse_regression$target <- factor(lumipulse_regression$target, levels = c(0, 1))

lumipulse_regression <- na.omit(lumipulse_regression)

numeric_vars <- c("AB40_lumipulse", "AB42_lumipulse", "Tau_lumipulse", "Age_onset_symptoms", "CMBS_deep", "CMBS_lobar")
categorical_vars <- setdiff(colnames(lumipulse_regression), c(numeric_vars, "target"))

# ---- Outlier Removal (Mahalanobis) ----
for(i in 1:2) {
  numeric_lum <- lumipulse_regression[, numeric_vars]
  M <- colMeans(numeric_lum)
  S <- cov(numeric_lum)
  d2 <- matrix(mahalanobis(numeric_lum, M, S))
  valid_indexes <- which(d2 <= 20)
  lumipulse_regression <- lumipulse_regression[valid_indexes, ]
}

# ---- Train / Test Split ----
set.seed(123)
train_idx <- createDataPartition(lumipulse_regression$target, p = 0.7, list = FALSE)
train_data <- lumipulse_regression[train_idx, ]
test_data  <- lumipulse_regression[-train_idx, ]

# ---- Data Augmentation (Applied ONLY to Train Data) ----
set.seed(123)
zeri_augmented <- augment_data(train_data, target_col = "target", 
                               numeric_cols = numeric_vars, categorical_cols = categorical_vars, 
                               n_new = 30)

train_data_augmented <- rbind(train_data, zeri_augmented)
train_data_augmented <- train_data_augmented[sample(nrow(train_data_augmented)), ]

# ---- Forward and Exhaustive Models (Trained on Augmented Train Data) ----
regfit.fwd <- regsubsets(target ~ ., data = train_data_augmented, nvmax = 8, method = "forward")
m1 <- glm(target ~ AB42_lumipulse + CMBS_deep + ICH + Antiplatelet + CMBS_lobar, 
          data = train_data_augmented, family = 'binomial')

regfit.full <- regsubsets(target ~ ., data = train_data_augmented, nvmax = 5)
m2 <- glm(target ~ AB42_lumipulse + CMBS_deep + ICH + Antiplatelet + Alchool, 
          data = train_data_augmented, family = 'binomial')
anova(m1, m2)

# ---- Lasso Regression ----
x_train <- model.matrix(target ~ ., data = train_data_augmented)[, -1]
y_train <- train_data_augmented$target
grid <- 10^seq(10, -2, length = 100)

lasso.mod <- glmnet(x_train, y_train, alpha = 1, lambda = grid, family = 'binomial')

set.seed(123)
cv.out <- cv.glmnet(x_train, y_train, alpha = 1, nfold = 10, lambda = grid, family = 'binomial') 
optlam.lasso <- cv.out$lambda.1se

coef.lasso <- predict(lasso.mod, s = optlam.lasso, type = 'coefficients')[1:20, ]

m3 <- glm(target ~ AB42_lumipulse + CMBS_deep + CMBS_lobar + Hypertension + ICH + Antiplatelet, 
          data = train_data_augmented, family = 'binomial')

# ---- Lasso with Interactions ----
x_train_int <- model.matrix(target ~ (AB42_lumipulse + CMBS_deep + CMBS_lobar + Hypertension + ICH + Antiplatelet)^2, 
                            data = train_data_augmented)[, -1]
lasso.mod_int <- glmnet(x_train_int, y_train, alpha = 1, lambda = grid, family = 'binomial')

set.seed(123)
cv.out_int <- cv.glmnet(x_train_int, y_train, alpha = 1, nfold = 10, lambda = grid, family = 'binomial') 
optlam.lasso_int <- cv.out_int$lambda.1se

coef.lasso_int <- predict(lasso.mod_int, s = optlam.lasso_int, type = 'coefficients')[1:20, ]

# Best Model
m4 <- glm(target ~ AB42_lumipulse + CMBS_lobar + ICH + Antiplatelet + CMBS_deep:Hypertension,
          data = train_data_augmented, family = 'binomial')
anova(m3, m4)

# ---- Dynamic Optimal Threshold via Youden's Index ----
prob_train <- predict(m4, type = "response")
roc_train <- roc(train_data_augmented$target, prob_train)
best_thresh <- coords(roc_train, "best", ret="threshold")$threshold[1]

# ---- Model Evaluation on Unseen TEST DATA ----
prob_test <- predict(m4, newdata = test_data, type = "response")
pred_class_test <- ifelse(prob_test >= best_thresh, 1, 0)

colors_true_test <- ifelse(test_data$target == 1, "red", "blue")
colors_pred_test <- ifelse(pred_class_test == 1, "darkgreen", "orange")

vars_to_plot <- c("AB42_lumipulse", "CMBS_lobar", "ICH", "Antiplatelet")

# Main effects plots
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (v in vars_to_plot) {
  plot(as.numeric(as.character(test_data[[v]])), prob_test, col = colors_true_test, pch = 19, 
       main = paste("Prob vs", v), xlab = v, ylab = "Predicted Probability")
}
par(mfrow = c(1, 1))

# Dedicated ggplot for the interaction term
ggplot(test_data, aes(x = factor(CMBS_deep), y = prob_test, fill = factor(Hypertension))) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(aes(color = factor(target)), position = position_dodge(0.75), size = 2) +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  labs(title = "Interaction Effect: CMBS_deep and Hypertension on Predicted Prob.",
       x = "CMBS_deep", y = "Predicted Probability", fill = "Hypertension", color = "True Target") +
  theme_minimal()

# ---- Evaluation Metrics (Test Set) ----
conf_obj_test <- confusionMatrix(factor(pred_class_test, levels = c(0, 1)), 
                                 factor(test_data$target, levels = c(0, 1)), 
                                 positive = '1')

conf_df_test <- as.data.frame(as.table(conf_obj_test$table))
colnames(conf_df_test) <- c("Prediction", "Reference", "Freq")
sensitivity_test <- round(conf_obj_test$byClass["Sensitivity"], 4)
specificity_test <- round(conf_obj_test$byClass["Specificity"], 4)

ggplot(conf_df_test, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), size = 6, color = "white") +
  scale_fill_gradient(low = "pink", high = "red") +
  labs(title = paste0("Test Set Confusion Matrix\nSensitivity: ", sensitivity_test, " | Specificity: ", specificity_test), 
       x = "Actual Target", y = "Predicted Class") +
  theme_minimal(base_size = 14)

roc_obj_test <- roc(test_data$target, prob_test)
plot(roc_obj_test, col = "red", main = "Test Set ROC Curve")
cat("Test AUC:", auc(roc_obj_test), "\n")

library(arm)
binnedplot(fitted(m4), rstandard(m4), main = "Binned Residual Plot (Train Set)", xlab = "Predicted Values", ylab = "Standardized Residuals")

# ---- Mixed-Effects Model (GLMM) ----
m_rand <- glmer(target ~ AB42_lumipulse + CMBS_lobar + ICH + Antiplatelet + CMBS_deep:Hypertension + (1 | Hypertension), 
                data = train_data_augmented, family = binomial)
summary(m_rand)

par(mfrow = c(2, 2))
plot(m_rand)
dotplot(ranef(m_rand))

sigma2_eps <- as.numeric(get_variance_residual(m_rand))
sigma2_b <- as.numeric(get_variance_random(m_rand))
PVRE <- (sigma2_b)^2 / (sigma2_b + sigma2_eps)^2

# ---- Plots for GLMM ----
plot_curve_continuous <- function(varname, data, model, levels_n = 100) {
  var_seq <- seq(min(data[[varname]], na.rm = TRUE), max(data[[varname]], na.rm = TRUE), length.out = levels_n)
  
  grid <- expand.grid(
    Hypertension = factor(c(0, 1), levels = levels(data$Hypertension)),
    AB42_lumipulse = median(data$AB42_lumipulse, na.rm = TRUE),
    CMBS_lobar = median(data$CMBS_lobar, na.rm = TRUE),
    CMBS_deep = median(data$CMBS_deep, na.rm = TRUE),
    ICH = factor(0, levels = levels(data$ICH)),
    Antiplatelet = factor(0, levels = levels(data$Antiplatelet))   )      grid <- grid[rep(1:nrow(grid), each = levels_n), ]   grid[[varname]] <- rep(var_seq, times = 2)   grid$`CMBS_deep:Hypertension` <- as.numeric(as.character(grid$Hypertension)) * grid$CMBS_deep
  grid$pred <- predict(model, newdata = grid, type = "response", allow.new.levels = TRUE)
  
  ggplot() +
    geom_jitter(data = data, aes_string(x = varname, y = "as.numeric(target)", color = "Hypertension"), height = 0.05, alpha = 0.4) +
    geom_line(data = grid, aes_string(x = varname, y = "pred", color = "Hypertension"), size = 1.2) +
    labs(x = varname, y = "Estimated probability of target = 1") +
    theme_minimal()
}

plot_curve_categorical <- function(varname, data, model) {
  grid <- expand.grid(
    Hypertension = factor(c(0, 1), levels = levels(data$Hypertension)),
    AB42_lumipulse = median(data$AB42_lumipulse, na.rm = TRUE),
    CMBS_lobar = median(data$CMBS_lobar, na.rm = TRUE),
    CMBS_deep = median(data$CMBS_deep, na.rm = TRUE),
    ICH = factor(0, levels = levels(data$ICH)),
    Antiplatelet = factor(0, levels = levels(data$Antiplatelet))   )      grid <- grid[rep(1:nrow(grid), each = length(levels(data[[varname]]))), ]   grid[[varname]] <- factor(rep(levels(data[[varname]]), times = nrow(grid) / length(levels(data[[varname]]))), levels = levels(data[[varname]]))   grid$`CMBS_deep:Hypertension` <- as.numeric(as.character(grid$Hypertension)) * grid$CMBS_deep
  grid$pred <- predict(model, newdata = grid, type = "response", allow.new.levels = TRUE)
  
  ggplot() +
    geom_jitter(data = data, aes_string(x = varname, y = "as.numeric(target)", color = "Hypertension"), width = 0.2, height = 0.05, alpha = 0.4) +
    geom_point(data = grid, aes_string(x = varname, y = "pred", color = "Hypertension"), size = 3, shape = 16, position = position_dodge(width = 0.4)) +
    labs(x = varname, y = "Estimated probability of target = 1") +
    theme_minimal()
}

plot_curve_continuous("AB42_lumipulse", train_data_augmented, m_rand)
plot_curve_continuous("CMBS_lobar", train_data_augmented, m_rand)
plot_curve_continuous("CMBS_deep", train_data_augmented, m_rand)
plot_curve_categorical("ICH", train_data_augmented, m_rand)
plot_curve_categorical("Antiplatelet", train_data_augmented, m_rand)
par(mfrow = c(1, 1))

# Evaluation on Test set for GLMM
prob_test_rand <- predict(m_rand, newdata = test_data, type = "response", allow.new.levels = TRUE)
pred_class_test_rand <- ifelse(prob_test_rand >= best_thresh, 1, 0)
cm_rand <- confusionMatrix(factor(pred_class_test_rand, levels = c(0,1)), factor(test_data$target, levels=c(0,1)), positive = '1')
print(cm_rand)

roc_obj_rand <- roc(test_data$target, prob_test_rand)
plot(roc_obj_rand, col = "blue", main = "Test ROC Curve - GLMM")
cat("GLMM Test AUC:", auc(roc_obj_rand), "\n")

# ---- GLMM Manual Cross Validation ----
set.seed(123)
folds <- createFolds(train_data_augmented$target, k = 5, list = TRUE, returnTrain = FALSE)

all_preds <- factor(levels = levels(train_data_augmented$target)) 
all_truth <- train_data_augmented$target

for (i in seq_along(folds)) {
  test_indices <- folds[[i]]
  cv_train_data <- train_data_augmented[-test_indices, ]
  cv_test_data <- train_data_augmented[test_indices, ]
  
  model <- glmer(target ~ AB42_lumipulse + CMBS_lobar + ICH + Antiplatelet + CMBS_deep:Hypertension + (1 | Hypertension),
                 data = cv_train_data, family = binomial)
  
  probs <- predict(model, newdata = cv_test_data, type = "response", allow.new.levels = TRUE)
  preds_class <- factor(ifelse(probs >= best_thresh, levels(cv_test_data$target)[2], levels(cv_test_data$target)[1]), levels = levels(cv_test_data$target))
  all_preds[test_indices] <- preds_class
  
  cat(sprintf("Confusion matrix fold %d:\n", i))
  print(confusionMatrix(preds_class, cv_test_data$target, positive = "1"))
  cat("\n-----------------------\n")
}

cat("Global confusion matrix on all folds (Train Data):\n")
print(confusionMatrix(all_preds, all_truth, positive = "1"))

# ---- Caret CV Evaluation (Performed on Train Data to avoid leakage) ----
train_data_augmented$target_factor <- factor(ifelse(train_data_augmented$target == 1, "Yes", "No"), levels = c("Yes", "No"))
formula <- target_factor ~ AB42_lumipulse + CMBS_lobar + ICH + Antiplatelet + CMBS_deep:Hypertension

train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

set.seed(123)
cv_model <- train(formula, data = train_data_augmented, method = "glm", family = "binomial", 
                  trControl = train_control, metric = "ROC")

cm_cv <- confusionMatrix(data = factor(cv_model$pred$pred, levels = c("No", "Yes")),
                         reference = factor(cv_model$pred$obs, levels = c("No", "Yes")), positive = "Yes")

conf_df_cv <- as.data.frame(as.table(cm_cv$table))
colnames(conf_df_cv) <- c("Prediction", "Reference", "Freq")

ggplot(conf_df_cv, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "black") + geom_text(aes(label = Freq), size = 6, color = "white") +
  scale_fill_gradient(low = "pink", high = "red") +
  labs(title = paste0("Confusion Matrix (5-Fold CV on Train)\nSensitivity: ", round(cm_cv$byClass["Sensitivity"], 4), " \vert{} Specificity: ", round(cm_cv$byClass["Specificity"], 4)), x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 14)

# ---- LOOCV Evaluation (Performed on Train Data) ----
train_control_loocv <- trainControl(method = "LOOCV", classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
set.seed(123)
loocv_model <- train(formula, data = train_data_augmented, method = "glm", family = "binomial", 
                     trControl = train_control_loocv, metric = "ROC")

cm_loocv <- confusionMatrix(data = factor(loocv_model$pred$pred, levels = c("No", "Yes")),
                            reference = factor(loocv_model$pred$obs, levels = c("No", "Yes")), positive = "Yes")

conf_df_loocv <- as.data.frame(as.table(cm_loocv$table))
colnames(conf_df_loocv) <- c("Prediction", "Reference", "Freq")

ggplot(conf_df_loocv, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "black") + geom_text(aes(label = Freq), size = 6, color = "white") +
  scale_fill_gradient(low = "pink", high = "red") +
  labs(title = paste0("Confusion Matrix (LOOCV on Train)\nSensitivity: ", round(cm_loocv$byClass["Sensitivity"], 4), " \vert{} Specificity: ", round(cm_loocv$byClass["Specificity"], 4)), x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 14)