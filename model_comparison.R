#!/usr/bin/env Rscript
# ===============================================
# 机器学习模型比较分析脚本
# 作者: 杜伟轩
# 功能:
#   1. 支持多种机器学习模型 (RF, XGBoost, SVM, Elastic Net)
#   2. 自动比较性能并生成可视化
#   3. 允许用户自定义并行核数
# ===============================================

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(ranger)
  library(xgboost)
  library(kernlab)
  library(ggplot2)
  library(patchwork)
  library(doParallel)
  library(parallel)
})

# ---------------- 主函数 ----------------
run_model_comparison <- function(data_file, out_dir = "model_comparison", n_cores = detectCores() - 1) {
  dir.create(out_dir, showWarnings = FALSE)

  # ---------------- 数据预处理 ----------------
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD, LSD, mC, Len), scale))
  X1 <- detagene2_scaled[, c("SSD", "LSD", "mC", "Len", "Sub")]
  y1 <- detagene2_scaled$Gene_mutation

  # Winsorize极端值
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(
    SSD = as.numeric(X1$SSD),
    LSD = as.numeric(X1$LSD),
    mC  = as.numeric(X1$mC),
    Len = as.numeric(X1$Len),
    Sub = factor(X1$Sub)
  )

  # 划分数据集
  train_idx <- createDataPartition(y_clean1, p = 0.8, list = FALSE)
  data_train <- data.frame(X_clean1[train_idx, ], Gene_mutation = y_clean1[train_idx])
  data_test  <- data.frame(X_clean1[-train_idx, ], Gene_mutation = y_clean1[-train_idx])
  data_train <- data_train %>% distinct(SSD, LSD, mC, Len, Sub, Gene_mutation, .keep_all = TRUE)

  # ---------------- 并行计算 ----------------
  cat("Using", n_cores, "CPU cores for parallel training...\n")
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    search = "grid",
    allowParallel = TRUE,
    savePredictions = "final"
  )

  set.seed(123)

  # ============= 模型训练部分 =============
  cat("Training Random Forest model...\n")
  rf_grid <- expand.grid(mtry = c(2, 3, 4),
                         splitrule = "variance",
                         min.node.size = c(5, 10, 15))
  rf_model <- train(
    Gene_mutation ~ ., data = data_train,
    method = "ranger", trControl = ctrl,
    tuneGrid = rf_grid, num.trees = 500,
    importance = "permutation", metric = "RMSE"
  )

  cat("Training XGBoost model...\n")
  xgb_grid <- expand.grid(
    nrounds = c(100, 200),
    max_depth = c(3, 6),
    eta = c(0.01, 0.1),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  )
  xgb_model <- train(
    Gene_mutation ~ ., data = data_train,
    method = "xgbTree", trControl = ctrl,
    tuneGrid = xgb_grid, metric = "RMSE",
    verbosity = 0
  )

  cat("Training SVM model...\n")
  svm_grid <- expand.grid(sigma = c(0.01, 0.1, 1),
                          C = c(0.1, 1, 10))
  svm_model <- train(
    Gene_mutation ~ ., data = data_train,
    method = "svmRadial", trControl = ctrl,
    tuneGrid = svm_grid, metric = "RMSE",
    preProcess = c("center", "scale")
  )

  cat("Training Elastic Net model...\n")
  enet_grid <- expand.grid(lambda = c(0.001, 0.01, 0.1),
                           fraction = seq(0.1, 0.9, length = 5))
  enet_model <- train(
    Gene_mutation ~ ., data = data_train,
    method = "enet", trControl = ctrl,
    tuneGrid = enet_grid, metric = "RMSE",
    preProcess = c("center", "scale")
  )

  stopCluster(cl)
  registerDoSEQ()

  # ---------------- 模型比较 ----------------
  models <- list(
    "Random Forest" = rf_model,
    "XGBoost"       = xgb_model,
    "SVM"           = svm_model,
    "Elastic Net"   = enet_model
  )

  test_predictions <- lapply(models, function(model) predict(model, newdata = data_test))

  test_metrics <- do.call(rbind, lapply(names(test_predictions), function(model_name) {
    pred <- test_predictions[[model_name]]
    obs  <- data_test$Gene_mutation
    data.frame(
      Model = model_name,
      RMSE = RMSE(pred, obs),
      Rsquared = R2(pred, obs),
      MAE = MAE(pred, obs)
    )
  }))

  write.csv(test_metrics, file.path(out_dir, "test_metrics_comparison.csv"), row.names = FALSE)
  saveRDS(models, file.path(out_dir, "trained_models.rds"))

  create_comparison_plots(models, test_metrics, out_dir, data_test)

  best_model <- test_metrics$Model[which.max(test_metrics$Rsquared)]
  cat("Best model:", best_model, "\n")
  cat("Test R-squared:", round(max(test_metrics$Rsquared), 4), "\n")

  return(list(models = models, test_metrics = test_metrics, best_model = best_model))
}

# ---------------- 可视化函数 ----------------
create_comparison_plots <- function(models, test_metrics, out_dir, data_test) {
  cat("Generating visualization results...\n")

  results <- resamples(models)
  results_df <- results$values

  # 1. Cross-validation R2
  p1 <- ggplot(results_df, aes(x = model, y = Rsquared)) +
    geom_boxplot(aes(fill = model), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "Cross-validation R2 Distribution", x = "Model", y = "R2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # 2. Cross-validation RMSE
  p2 <- ggplot(results_df, aes(x = model, y = RMSE)) +
    geom_boxplot(aes(fill = model), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "Cross-validation RMSE Distribution", x = "Model", y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # 3. Test set R2
  p3 <- ggplot(test_metrics, aes(x = reorder(Model, Rsquared), y = Rsquared, fill = Model)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = round(Rsquared, 3)), vjust = -0.5) +
    labs(title = "Test Set R2 Comparison", x = "Model", y = "R2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  model_comparison_plot <- (p1 + p2 + p3) + plot_layout(ncol = 3)
  ggsave(file.path(out_dir, "model_comparison.png"), model_comparison_plot,
         width = 15, height = 5, dpi = 300)

  # Save summary table
  performance_table <- test_metrics %>%
    arrange(desc(Rsquared)) %>%
    mutate(across(where(is.numeric), round, 4))
  write.csv(performance_table, file.path(out_dir, "model_performance_summary.csv"), row.names = FALSE)

  # Predicted vs Observed
  predictions_combined <- data.frame()
  for (model_name in names(models)) {
    pred <- predict(models[[model_name]], data_test)
    temp_df <- data.frame(
      Model = model_name,
      Observed = data_test$Gene_mutation,
      Predicted = pred
    )
    predictions_combined <- rbind(predictions_combined, temp_df)
  }

  p_scatter <- ggplot(predictions_combined, aes(x = Observed, y = Predicted, color = Model)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Model, ncol = 2) +
    labs(title = "Predicted vs Observed", x = "Observed", y = "Predicted") +
    theme_minimal()

  ggsave(file.path(out_dir, "prediction_scatter.png"), p_scatter, width = 10, height = 8, dpi = 300)

  cat("Visualization results saved to:", out_dir, "\n")
}

# ---------------- 脚本入口 ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript model_comparison.R <data_file.csv> [output_dir] [n_cores]\n")
} else {
  data_file <- args[1]
  out_dir <- ifelse(length(args) > 1, args[2], "model_comparison")
  n_cores <- ifelse(length(args) > 2, as.numeric(args[3]), detectCores() - 1)
  run_model_comparison(data_file, out_dir, n_cores)
}
