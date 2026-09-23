# =====================================================================
# ENSEMBLE MODELS 
# =====================================================================

library(MASS)
library(caret)
library(recipes)
library(dplyr)
library(pROC)
library(caretEnsemble)
library(VIM)
library(GGally)
library(e1071)
library(boot)
library(ggplot2)
library(reshape2)
library(glmnet)
library(car)

# ---- Data Import & Preprocessing ----
load('data2.RData')

data.ens <- data.frame(
  target = data2$`Radiological definition of the pathology (CAA / CAA-DPA / DPA)_RM_REC`,
  AB40_plasma = data2$`Aβ40(pg/ml)PLASMA`,
  AB42_plasma = data2$`Aβ42 (pg/ml)PLASMA`,
  Tau_plasma = data2$`Tau (pg/ml)PLASMA`,
  CMB_lobar = data2$`CMBs (lobar) _RM_REC`,
  CMB_deep = data2$`CMBs (deep)_RM_REC`,
  Dyslipidemia_RISK = data2$Dyslipidemia_RISK,   ICH = data2$`ICH (1/0)_RECRUIT`, 
  Antiplatelet = data2$`antiplatelet (type)_TERAPIA_REC`
)

data.ens$target <- ifelse(data.ens$target %in% c("CAA", "CAA-DPA", "iCAA"), 1, 0)
data.ens$target <- factor(data.ens$target, levels = c(0, 1), labels = c("no", "yes"))

data.ens <- na.omit(data.ens)

# ---- Outlier Removal (Mahalanobis) ----
numeric_vars_raw <- c("AB40_plasma", "AB42_plasma", "Tau_plasma", "CMB_deep", "CMB_lobar")
numeric_ens <- data.ens[, numeric_vars_raw]

M <- colMeans(numeric_ens)
S <- cov(numeric_ens)
d2 <- matrix(mahalanobis(numeric_ens, M, S))

data.ens <- data.ens[which(d2 <= 20), ]

# ---- Feature Transformation ----
data.ens$CMB_deep_log <- log(1 + data.ens$CMB_deep)

numeric_vars <- c("AB40_plasma", "AB42_plasma", "Tau_plasma", "CMB_deep_log", "CMB_lobar")
categorical_vars <- c("Dyslipidemia_RISK", "ICH", "Antiplatelet")

# ---- Train/Test Split ----
set.seed(123)
index_train <- createDataPartition(data.ens$target, p = 0.9, list = FALSE)
data_train <- data.ens[index_train, ]
data_test  <- data.ens[-index_train, ]

# ---- Data Augmentation ----
augment_class0 <- function(data, target_col = "target", numeric_cols = NULL, categorical_cols = NULL, n_new = 20, k_range = c(2, 5), noise_sd = 0.05, k_nn = 3) {
  class0 <- subset(data, data[[target_col]] == "no")
  n_class0 <- nrow(class0)
  
  if (n_class0 < max(k_range)) stop("Class 0 has fewer observations than the maximum value of k.")
  
  augmented <- data.frame()
  
  for (i in 1:n_new) {
    k <- sample(k_range[1]:k_range[2], 1)
    idx <- sample(1:n_class0, k, replace = FALSE)
    selected_numeric <- class0[idx, numeric_cols, drop = FALSE]
    
    weights <- runif(k) / sum(runif(k))
    combined_numeric <- colSums(t(weights * t(as.matrix(selected_numeric))))
    noisy_combined <- combined_numeric + rnorm(length(combined_numeric), mean = 0, sd = noise_sd)
    
    distances <- apply(class0[, numeric_cols, drop = FALSE], 1, function(row) sqrt(sum((row - noisy_combined)^2)))
    nearest_cat <- class0[order(distances)[1:k_nn], categorical_cols, drop = FALSE]
    assigned_cat <- sapply(nearest_cat, function(col) names(which.max(table(col))), simplify = TRUE)
    
    new_point_df <- as.data.frame(as.list(c(noisy_combined, assigned_cat)), stringsAsFactors = FALSE)
    new_point_df[[target_col]] <- "no"
    
    for (v in numeric_cols) new_point_df[[v]] <- as.numeric(new_point_df[[v]])
    for (v in categorical_cols) new_point_df[[v]] <- factor(new_point_df[[v]], levels = levels(data[[v]]))
    
    augmented <- rbind(augmented, new_point_df)
  }
  
  return(augmented[, colnames(data)])
}

set.seed(123)
zeri_augmented <- augment_class0(data_train, numeric_cols = numeric_vars, categorical_cols = categorical_vars, n_new = 30)
data_train_aug <- rbind(data_train, zeri_augmented)
data_train_aug <- data_train_aug[sample(nrow(data_train_aug)), ]

# ---- Scaling ----
preproc <- preProcess(data_train_aug[, numeric_vars], method = c("center", "scale"))
data_train_aug[, numeric_vars] <- predict(preproc, data_train_aug[, numeric_vars])
data_test[, numeric_vars] <- predict(preproc, data_test[, numeric_vars])

# Remove the raw CMB_deep to avoid collinearity with CMB_deep_log during modeling
lasso_train_data <- data_train_aug %>% dplyr::select(-CMB_deep)
lasso_test_data <- data_test %>% dplyr::select(-CMB_deep)

# ---- Bootstrap Variable Selection (Lasso Stability) ----
set.seed(123)
B <- 100 

# Create design matrix with all two-way interactions
X_train_raw <- model.matrix(target ~ .^2 - 1, data = lasso_train_data)
# Make names syntactically valid for R formulas (replaces ':' with '.')
valid_names <- make.names(colnames(X_train_raw))
colnames(X_train_raw) <- valid_names
y_train <- lasso_train_data$target

selected_count <- setNames(rep(0, ncol(X_train_raw)), colnames(X_train_raw))

for (b in 1:B) {
  idx <- sample(1:nrow(X_train_raw), replace = TRUE)
  cvfit <- cv.glmnet(X_train_raw[idx, ], y_train[idx], family = "binomial", alpha = 1, standardize = FALSE, type.measure = "class")
  
  coef_b <- coef(cvfit, s = "lambda.min")
  vars_b <- rownames(coef_b)[coef_b[, 1] != 0][-1]
  selected_count[vars_b] <- selected_count[vars_b] + 1
}

selection_freq_df <- data.frame(Variable = names(selected_count), SelectionRate = (selected_count / B) * 100) %>% 
  arrange(desc(SelectionRate))

ggplot(selection_freq_df[1:20, ], aes(x = reorder(Variable, SelectionRate), y = SelectionRate, fill = SelectionRate)) +
  geom_col() + coord_flip() + scale_fill_gradient(low = "#FFCCCC", high = "#990000") +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray", linewidth = 0.8) +
  labs(title = "Variable Selection Frequencies (Top 20 Interactions)", x = "", y = "", fill = "Selection Rate") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")

# ---- Final Model Data Preparation ----
selected_vars <- selection_freq_df %>% filter(SelectionRate > 80) %>% pull(Variable)

X_train_df <- as.data.frame(X_train_raw)
X_test_raw <- model.matrix(target ~ .^2 - 1, data = lasso_test_data)
colnames(X_test_raw) <- make.names(colnames(X_test_raw))
X_test_df <- as.data.frame(X_test_raw)

data_train_final <- X_train_df[, selected_vars, drop = FALSE]
data_train_final$target <- y_train

data_test_final <- X_test_df[, selected_vars, drop = FALSE]
data_test_final$target <- lasso_test_data$target

# ---- Repeated CV (50x20) on Selected Variables ----
set.seed(123)
n_repeats <- 50
n_folds   <- 20

accuracy_all <- c(); sensitivity_all <- c(); specificity_all <- c(); f1_all <- c()

for (rep in 1:n_repeats) {
  yes_data <- data_train_final %>% filter(target == "yes") %>% slice_sample(prop = 1)
  no_data  <- data_train_final %>% filter(target == "no")  %>% slice_sample(prop = 1)
  
  folds_yes <- cut(seq_len(nrow(yes_data)), breaks = n_folds, labels = FALSE)
  folds_no  <- cut(seq_len(nrow(no_data)),  breaks = n_folds, labels = FALSE)
  
  for (i in 1:n_folds) {
    fold_test <- bind_rows(yes_data[which(folds_yes == i), ], no_data[which(folds_no == i), ]) %>% arrange(row_number())
    fold_train <- bind_rows(yes_data[which(folds_yes != i), ], no_data[which(folds_no != i), ]) %>% arrange(row_number())
    
    model <- suppressWarnings(glm(target ~ ., data = fold_train, family = "binomial"))
    
    preds <- ifelse(predict(model, newdata = fold_test, type = "response") > 0.5, "yes", "no") %>% factor(levels = c("no", "yes"))
    cm <- confusionMatrix(preds, fold_test$target, positive = "yes")
    
    precision <- cm$byClass["Precision"]
    recall    <- cm$byClass["Sensitivity"]
    
    accuracy_all <- c(accuracy_all, cm$overall["Accuracy"])
    sensitivity_all <- c(sensitivity_all, cm$byClass["Sensitivity"])
    specificity_all <- c(specificity_all, cm$byClass["Specificity"])
    f1_all <- c(f1_all, ifelse((precision + recall) == 0, NA, 2 * precision * recall / (precision + recall)))
  }
}

cv.results <- function(acc, sens, spec, f1, f, r) {
  cat("╔══════════════════════════════════════════════════════╗\n")
  cat("║               Cross-Validation Results               ║\n")
  cat("╚══════════════════════════════════════════════════════╝\n")
  cat(sprintf("Accuracy   : %.4f ± %.4f\nSensitivity: %.4f ± %.4f\nSpecificity: %.4f ± %.4f\nF1 Score   : %.4f ± %.4f\n", 
              mean(acc, na.rm = TRUE), sd(acc, na.rm = TRUE), mean(sens, na.rm = TRUE), sd(sens, na.rm = TRUE),
              mean(spec, na.rm = TRUE), sd(spec, na.rm = TRUE), mean(f1, na.rm = TRUE), sd(f1, na.rm = TRUE)))
}
cv.results(accuracy_all, sensitivity_all, specificity_all, f1_all, n_folds, n_repeats)

# ---- LOOCV on Selected Variables ----
set.seed(123)
accuracy_loocv <- c(); sensitivity_loocv <- c(); specificity_loocv <- c(); f1_loocv <- c()

for (i in 1:nrow(data_train_final)) {
  fold_test <- data_train_final[i, , drop = FALSE]
  fold_train <- data_train_final[-i, , drop = FALSE]
  
  model <- suppressWarnings(glm(target ~ ., data = fold_train, family = "binomial"))
  
  preds <- ifelse(predict(model, newdata = fold_test, type = "response") > 0.5, "yes", "no") %>% factor(levels = c("no", "yes"))
  cm <- confusionMatrix(preds, fold_test$target, positive = "yes")
  
  precision <- cm$byClass["Precision"]
  recall    <- cm$byClass["Sensitivity"]
  
  accuracy_loocv <- c(accuracy_loocv, cm$overall["Accuracy"])
  sensitivity_loocv <- c(sensitivity_loocv, cm$byClass["Sensitivity"])
  specificity_loocv <- c(specificity_loocv, cm$byClass["Specificity"])
  f1_loocv <- c(f1_loocv, ifelse((precision + recall) == 0, NA, 2 * precision * recall / (precision + recall)))
}

loocv.results <- function(acc, sens, spec, f1) {
  cat("╔═══════════════════════════════════════════════════╗\n")
  cat("║                   LOOCV Results                   ║\n")
  cat("╚═══════════════════════════════════════════════════╝\n")
  cat(sprintf("Accuracy   : %.4f ± %.4f\nSensitivity: %.4f ± %.4f\nSpecificity: %.4f ± %.4f\nF1 Score   : %.4f ± %.4f\n", 
              mean(acc, na.rm = TRUE), sd(acc, na.rm = TRUE), mean(sens, na.rm = TRUE), sd(sens, na.rm = TRUE),
              mean(spec, na.rm = TRUE), sd(spec, na.rm = TRUE), mean(f1, na.rm = TRUE), sd(f1, na.rm = TRUE)))
}
loocv.results(accuracy_loocv, sensitivity_loocv, specificity_loocv, f1_loocv)

# ---- Plot 95% CI for LOOCV metrics ----
df_ci <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity"),
  Mean   = c(mean(accuracy_loocv, na.rm = TRUE), mean(sensitivity_loocv, na.rm = TRUE), mean(specificity_loocv, na.rm = TRUE)),
  SE     = c(sd(accuracy_loocv, na.rm = TRUE), sd(sensitivity_loocv, na.rm = TRUE), sd(specificity_loocv, na.rm = TRUE)) / sqrt(length(accuracy_loocv))
)
df_ci$Lower <- df_ci$Mean - 1.96 * df_ci$SE
df_ci$Upper <- pmin(df_ci$Mean + 1.96 * df_ci$SE, 1)

ggplot(df_ci, aes(x = Metric, y = Mean)) +
  geom_point(size = 3, color = "red4") + geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, lwd = 1, color = "red") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray") + ylim(0, 1) +
  theme_minimal(base_size = 13) + labs(title = "95% CI for LOOCV Metrics", y = "", x = "")

# ---- Beta Coefficients Confidence Intervals ----
model_final <- glm(target ~ ., data = data_train_final, family = "binomial")
ci <- confint.default(model_final, level = 0.95)

df_plot <- data.frame(Term = factor(names(coef(model_final)), levels = names(coef(model_final))[order(coef(model_final))]),
                      Beta = coef(model_final), CI_Lower = ci[, 1], CI_Upper = ci[, 2])

ggplot(df_plot, aes(x = Beta, y = Term)) +
  geom_point(color = "red4", size = 2) + geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, lwd = 1, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + theme_minimal(base_size = 13) +
  labs(title = "CI for Betas", x = "Betas", y = "Predictors") + theme(panel.grid.minor = element_blank())

# ---- Bootstrap Beta Distribution ----
set.seed(123)
coeff_matrix <- matrix(NA, nrow = 200, ncol = length(coef(model_final)), dimnames = list(NULL, names(coef(model_final))))

for (b in 1:200) {
  boot_data <- data_train_final[sample(1:nrow(data_train_final), replace = TRUE), ]
  model_boot <- tryCatch(glm(target ~ ., data = boot_data, family = "binomial"), error = function(e) NULL)
  if (!is.null(model_boot)) coeff_matrix[b, ] <- coef(model_boot)
}

coeff_df <- subset(melt(as.data.frame(coeff_matrix[complete.cases(coeff_matrix), ]), variable.name = "Coefficient", value.name = "Beta"), abs(Beta) <= 1e3)

ggplot(coeff_df, aes(x = Coefficient, y = Beta)) +
  geom_boxplot(fill = "#69b3a2", alpha = 0.7) + geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() + labs(title = "Bootstrap Beta Coefficients Distribution", x = "Coefficient", y = "Beta Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Ensemble Learning (Bagging) ----
set.seed(123)
n_bootstraps <- 100
models_list <- vector("list", n_bootstraps)
predictions_list <- vector("list", n_bootstraps)

for (i in seq_len(n_bootstraps)) {
  boot_data <- data_train_final[sample(seq_len(nrow(data_train_final)), replace = TRUE), ]
  model <- suppressWarnings(glm(target ~ ., data = boot_data, family = "binomial"))
  models_list[[i]] <- model
  predictions_list[[i]] <- predict(model, newdata = data_test_final, type = "response")
}

# ---- Bagging Evaluation & Metric CIs ----
avg_probs <- rowMeans(do.call(cbind, predictions_list))

stats_fun <- function(data, indices) {
  pred_bs <- factor(ifelse(avg_probs[indices] > 0.5, "yes", "no"), levels = c("no", "yes"))
  cm_bs <- confusionMatrix(pred_bs, data_test_final$target[indices], positive = "yes")
  
  precision_bs <- cm_bs$byClass["Precision"]
  recall_bs    <- cm_bs$byClass["Sensitivity"]
  f1_bs <- ifelse(sum(c(precision_bs, recall_bs), na.rm = TRUE) == 0, 0, 2 * prod(c(precision_bs, recall_bs)) / sum(c(precision_bs, recall_bs)))
  
  c(AUC = auc(roc(data_test_final$target[indices], avg_probs[indices])), Accuracy = cm_bs$overall["Accuracy"],
    Sensitivity = recall_bs, Specificity = cm_bs$byClass["Specificity"], F1 = f1_bs)
}

set.seed(123)
boot_out <- boot(data = seq_along(avg_probs), statistic = stats_fun, R = 2000)

cat("╔═══════════════════════════════════════════════╗\n║       Bootstrap Confidence Intervals (95%)    ║\n╚═══════════════════════════════════════════════╝\n")
metriche <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1")
for (i in seq_along(metriche)) {
  ci <- boot.ci(boot_out, type = "perc", index = i)
  cat(sprintf("%-12s: %.3f (%.3f – %.3f)\n", metriche[i], boot_out$t0[i], ci$percent[4], ci$percent[5]))
}

boot_long <- melt(setNames(as.data.frame(boot_out$t), metriche), variable.name = "Metric", value.name = "Value")
ggplot(boot_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.7, color = "black") + theme_minimal(base_size = 14) + theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2") + geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(title = "Bootstrap Performance Metrics (Bagged Model Test Set)", y = "Score", x = NULL)

# ---- Bagged Predictions Variance Analysis ----
probs_var <- apply(do.call(cbind, predictions_list), 1, var)

ggplot(data.frame(prob_var = probs_var), aes(x = prob_var)) +
  geom_histogram(binwidth = 0.0005, fill = "steelblue", color = "white") +
  labs(title = "Variance distribution among bagged predictions", x = "Variance", y = "Frequency")

# ---- Analyze Coefficient Variability in Bagged Models ----
coef_long <- melt(data.frame(t(sapply(models_list, coef)), bootstrap = 1:100), id.vars = "bootstrap", variable.name = "term", value.name = "estimate")

ggplot(coef_long, aes(x = term, y = estimate)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coefficient variability across bagged models", x = "Predictors", y = "Estimate")