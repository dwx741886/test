#!/usr/bin/env Rscript
## 模型比较脚本
suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(ranger)
  library(xgboost)
  library(kernlab)
  library(ggplot2)
  library(patchwork)
  library(doParallel)
})

run_model_comparison <- function(data_file, out_dir = "model_comparison") {
  dir.create(out_dir, showWarnings = FALSE)
  
  # 数据预处理（与主脚本一致）
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD, LSD, mC, Len), scale))
  X1 <- detagene2_scaled[, c("SSD", "LSD", "mC", "Len", "Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(
    SSD = as.numeric(X1$SSD),
    LSD = as.numeric(X1$LSD),
    mC = as.numeric(X1$mC),
    Len = as.numeric(X1$Len),
    Sub = factor(X1$Sub)
  )
  
  # 数据集划分
  train_idx <- createDataPartition(y_clean1, p = 0.8, list = FALSE)
  data_train <- data.frame(X_clean1[train_idx, ], Gene_mutation = y_clean1[train_idx])
  data_test <- data.frame(X_clean1[-train_idx, ], Gene_mutation = y_clean1[-train_idx])
  data_train <- data_train %>% distinct(SSD, LSD, mC, Len, Sub, Gene_mutation, .keep_all = TRUE)
  
  # 设置并行计算
  cl <- makeCluster(7)
  registerDoParallel(cl)
  
  # 统一的训练控制参数
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    search = "grid",
    allowParallel = TRUE,
    savePredictions = "final"
  )
  
  # 1. 随机森林 (ranger)
  cat("训练随机森林模型...\n")
  rf_grid <- expand.grid(
    mtry = c(2, 3, 4),
    splitrule = "variance",
    min.node.size = c(5, 10, 15)
  )
  
  set.seed(123)
  rf_model <- train(
    Gene_mutation ~ .,
    data = data_train,
    method = "ranger",
    trControl = ctrl,
    tuneGrid = rf_grid,
    num.trees = 500,
    importance = "permutation",
    metric = "RMSE"
  )
  
  # 2. 梯度提升树 (XGBoost)
  cat("训练XGBoost模型...\n")
  xgb_grid <- expand.grid(
    nrounds = c(100, 200),
    max_depth = c(3, 6),
    eta = c(0.01, 0.1),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  )
  
  set.seed(123)
  xgb_model <- train(
    Gene_mutation ~ .,
    data = data_train,
    method = "xgbTree",
    trControl = ctrl,
    tuneGrid = xgb_grid,
    metric = "RMSE",
    verbosity = 0
  )
  
  # 3. 支持向量机 (SVM)
  cat("训练SVM模型...\n")
  svm_grid <- expand.grid(
    sigma = c(0.01, 0.1, 1),
    C = c(0.1, 1, 10)
  )
  
  set.seed(123)
  svm_model <- train(
    Gene_mutation ~ .,
    data = data_train,
    method = "svmRadial",
    trControl = ctrl,
    tuneGrid = svm_grid,
    metric = "RMSE",
    preProcess = c("center", "scale")
  )
  
  # 4. 弹性网络回归 (作为线性方法基准)
  cat("训练弹性网络模型...\n")
  enet_grid <- expand.grid(
    lambda = c(0.001, 0.01, 0.1),
    fraction = seq(0.1, 0.9, length = 5)
  )
  
  set.seed(123)
  enet_model <- train(
    Gene_mutation ~ .,
    data = data_train,
    method = "enet",
    trControl = ctrl,
    tuneGrid = enet_grid,
    metric = "RMSE",
    preProcess = c("center", "scale")
  )
  
  stopCluster(cl)
  registerDoSEQ()
  
  # 模型性能比较
  models <- list(
    "Random Forest" = rf_model,
    "XGBoost" = xgb_model,
    "SVM" = svm_model,
    "Elastic Net" = enet_model
  )
  
  # 提取交叉验证结果
  results <- resamples(models)
  
  # 测试集预测
  test_predictions <- lapply(models, function(model) {
    predict(model, newdata = data_test)
  })
  
  test_metrics <- do.call(rbind, lapply(names(test_predictions), function(model_name) {
    pred <- test_predictions[[model_name]]
    obs <- data_test$Gene_mutation
    data.frame(
      Model = model_name,
      RMSE = RMSE(pred, obs),
      Rsquared = R2(pred, obs),
      MAE = MAE(pred, obs)
    )
  }))
  
  # 保存结果
  write.csv(test_metrics, file.path(out_dir, "test_metrics_comparison.csv"), row.names = FALSE)
  saveRDS(models, file.path(out_dir, "trained_models.rds"))
  
  # 可视化比较
  create_comparison_plots(results, test_metrics, out_dir)
  
  # 返回最佳模型
  best_model <- names(which.min(test_metrics$RMSE))
  cat("最佳模型:", best_model, "\n")
  cat("测试集R²:", round(test_metrics$Rsquared[test_metrics$Model == best_model], 4), "\n")
  
  return(list(models = models, test_metrics = test_metrics, best_model = best_model))
}

# 创建比较图
create_comparison_plots <- function(results, test_metrics, out_dir) {
  
  # 1. 交叉验证性能比较
  p1 <- ggplot(results, aes(x = model, y = Rsquared)) +
    geom_boxplot(aes(fill = model), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "交叉验证R²分布", x = "模型", y = "R²") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(results, aes(x = model, y = RMSE)) +
    geom_boxplot(aes(fill = model), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "交叉验证RMSE分布", x = "模型", y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 2. 测试集性能比较
  p3 <- ggplot(test_metrics, aes(x = reorder(Model, Rsquared), y = Rsquared, fill = Model)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = round(Rsquared, 3)), vjust = -0.5) +
    labs(title = "测试集R²比较", x = "模型", y = "R²") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 3. 预测值 vs 观测值散点图
  model_comparison_plot <- p1 + p2 + p3 +
    plot_layout(ncol = 3) &
    theme(legend.position = "none")
  
  ggsave(file.path(out_dir, "model_comparison.png"), model_comparison_plot, 
         width = 15, height = 5, dpi = 300)
  
  # 4. 性能指标汇总表
  performance_table <- test_metrics %>%
    arrange(desc(Rsquared)) %>%
    mutate(across(where(is.numeric), round, 4))
  
  write.csv(performance_table, file.path(out_dir, "model_performance_summary.csv"), row.names = FALSE)
  
  cat("模型比较完成！结果保存至:", out_dir, "\n")
}

# 执行脚本
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("用法: Rscript model_comparison.R <data_file.csv> [output_dir]\n")
} else {
  data_file <- args[1]
  out_dir <- ifelse(length(args) > 1, args[2], "model_comparison")
  run_model_comparison(data_file, out_dir)
}
