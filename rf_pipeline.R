#!/usr/bin/env Rscript
# ===============================================
# 机器学习验证分析脚本
# 作者: 杜伟轩
# 功能:
#   1. 支持数据清洗、特征筛选
#   2. 支持随机森林模型
#   2. 10-fold、验证集验证
# ===============================================
##使用Rscript rf_pipeline.R ALLdistanceHC9.csv results_out/
##  1. 加载必要包 
suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(ranger)
  library(Boruta)
  library(pdp)
  library(ggplot2)
  library(doParallel)
  library(tidyverse)
  library(car)
})

## 2. main function
run_rf_pipeline <- function(data_file, out_dir = "results_out") {
  dir.create(out_dir, showWarnings = FALSE)
  
  # 数据读取与预处理 
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD, LSD, mC, Len), scale))
  X1 <- detagene2_scaled[, c("SSD", "LSD", "mC", "Len", "Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  
  # Winsorize
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(
    SSD = as.numeric(X1$SSD),
    LSD = as.numeric(X1$LSD),
    mC = as.numeric(X1$mC),
    Len = as.numeric(X1$Len),
    Sub = factor(X1$Sub)
  )
  
  # 分层抽样
  train_idx <- createDataPartition(y_clean1, p = 0.8, list = FALSE)
  data_train1 <- data.frame(X_clean1[train_idx, ], Gene_mutation = y_clean1[train_idx])
  data_ver1   <- data.frame(X_clean1[-train_idx, ], Gene_mutation = y_clean1[-train_idx])
  
  # 去重
  data_clean <- data_train1 %>% distinct(SSD, LSD, mC, Len, Sub, Gene_mutation, .keep_all = TRUE)
  cat("训练集样本数:", nrow(data_clean), "\n")
  
  # 特征选择 (Boruta) 
  cat("运行 Boruta 特征选择...\n")
  boruta <- Boruta(x = data_clean[, c("SSD","LSD","mC","Len","Sub")],
                   y = data_clean$Gene_mutation, maxRuns = 500, pValue = 0.01, mcAdj = TRUE)
  saveRDS(boruta, file.path(out_dir, "boruta_results.rds"))
  
  # 共线性分析 
  vif_values <- vif(lm(Gene_mutation ~ ., data = data_clean))
  write.csv(vif_values, file.path(out_dir, "vif_values.csv"))
  
  # 随机森林建模
  cat("开始超参调优...\n")
  cl <- makeCluster(7); registerDoParallel(cl)
  tuneGrid <- expand.grid(mtry = c(2, 3, 4), splitrule = "variance", min.node.size = c(5, 10, 15))
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "grid", allowParallel = TRUE)
  
  set.seed(123)
  rf_tune <- train(Gene_mutation ~ ., data = data_clean, method = "ranger",
                   trControl = ctrl, tuneGrid = tuneGrid, num.trees = 500,
                   importance = "permutation", metric = "RMSE")
  stopCluster(cl); registerDoSEQ()
  saveRDS(rf_tune, file.path(out_dir, "rf_tune.rds"))
  
  # 最终模型 
  rf_model <- ranger(Gene_mutation ~ ., data = data_clean, num.trees = 500,
                     mtry = 3, min.node.size = 10, importance = "permutation",
                     oob.error = TRUE, keep.inbag = TRUE, write.forest = TRUE)
  saveRDS(rf_model, file.path(out_dir, "rf_model.rds"))
  
  #  十折交叉验证
  ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = "final")
  cv_model <- train(Gene_mutation ~ ., data = data_clean, method = "ranger",
                    num.trees = 500, importance = "permutation",
                    trControl = ctrl, tuneGrid = data.frame(mtry=3,splitrule="variance",min.node.size=10))
  
  cv_r2_mean <- mean(cv_model$resample$Rsquared)
  cv_r2_sd   <- sd(cv_model$resample$Rsquared)
  cat("十折交叉验证 R²:", round(cv_r2_mean,3), "±", round(cv_r2_sd,3), "\n")
  
  # 测试集验证
  pred_values <- predict(rf_model, data = data_ver1)$predictions
  test_results <- postResample(pred = pred_values, obs = data_ver1$Gene_mutation)
  cat("测试集结果:\n"); print(test_results)
  
  # 绘图 
  p <- ggplot(cv_model$resample, aes(x = Rsquared)) +
    geom_histogram(fill = "steelblue", bins = 10) +
    geom_vline(xintercept = test_results[["Rsquared"]], color = "red", linetype = "dashed") +
    labs(title = "R² distribution (10-fold CV)", x = "R²", y = "Counts",
         subtitle = paste0("Mean = ", round(cv_r2_mean, 3), " (SD=", round(cv_r2_sd,3), ")"),
         caption = paste("Validation R² =", round(test_results[["Rsquared"]],3))) +
    theme_minimal()
  ggsave(file.path(out_dir, "cv_r2_distribution.png"), p, width=6, height=4)
}

# 3. 运行脚本 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript rf_pipeline.R <data_file.csv> <output_dir>\n")
} else {
  run_rf_pipeline(args[1], args[2])
}
