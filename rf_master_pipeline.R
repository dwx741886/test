#!/usr/bin/env Rscript
## 用法：
##  Rscript rf_master_pipeline.R <mode> <data_file.csv> <output_dir> [额外参数]
##  mode 可选：rf_pipeline / rf_shap / rf_interaction / permutation_test
##  permutation_test 模式需要额外参数: <n_perm> <grid_points> <num_cores>

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
  library(iml)
  library(fastshap)
  library(shapviz)
  library(ggbeeswarm)
  library(patchwork)
  library(foreach)
})

##1. rf_pipeline (建模 + 特征选择 + CV) 
run_rf_pipeline <- function(data_file, out_dir) {
  dir.create(out_dir, showWarnings = FALSE)
  
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD, LSD, mC, Len), scale))
  X1 <- detagene2_scaled[, c("SSD", "LSD", "mC", "Len", "Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(
    SSD = as.numeric(X1$SSD), LSD = as.numeric(X1$LSD),
    mC = as.numeric(X1$mC), Len = as.numeric(X1$Len),
    Sub = factor(X1$Sub)
  )
  
  train_idx <- createDataPartition(y_clean1, p = 0.8, list = FALSE)
  data_train1 <- data.frame(X_clean1[train_idx, ], Gene_mutation = y_clean1[train_idx])
  data_ver1   <- data.frame(X_clean1[-train_idx, ], Gene_mutation = y_clean1[-train_idx])
  
  data_clean <- data_train1 %>% distinct(SSD, LSD, mC, Len, Sub, Gene_mutation, .keep_all = TRUE)
  cat("训练集样本数:", nrow(data_clean), "\n")
  
  # Boruta
  cat("运行 Boruta 特征选择...\n")
  boruta <- Boruta(x = data_clean[, c("SSD","LSD","mC","Len","Sub")],
                   y = data_clean$Gene_mutation, maxRuns = 500, pValue = 0.01, mcAdj = TRUE)
  saveRDS(boruta, file.path(out_dir, "boruta_results.rds"))
  
  # VIF
  vif_values <- vif(lm(Gene_mutation ~ ., data = data_clean))
  write.csv(vif_values, file.path(out_dir, "vif_values.csv"))
  
  # 调参
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
  
  rf_model <- ranger(Gene_mutation ~ ., data = data_clean, num.trees = 500,
                     mtry = 3, min.node.size = 10, importance = "permutation",
                     oob.error = TRUE, keep.inbag = TRUE, write.forest = TRUE)
  saveRDS(rf_model, file.path(out_dir, "rf_model.rds"))
  
  # CV
  ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = "final")
  cv_model <- train(Gene_mutation ~ ., data = data_clean, method = "ranger",
                    num.trees = 500, importance = "permutation",
                    trControl = ctrl, tuneGrid = data.frame(mtry=3,splitrule="variance",min.node.size=10))
  
  cv_r2_mean <- mean(cv_model$resample$Rsquared)
  cv_r2_sd   <- sd(cv_model$resample$Rsquared)
  cat("十折交叉验证 R²:", round(cv_r2_mean,3), "±", round(cv_r2_sd,3), "\n")
  
  # 测试集
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

##2. rf_shap_pipeline (SHAP 分析) 
run_rf_shap <- function(data_file, out_dir) {
  dir.create(out_dir, showWarnings = FALSE)
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD,LSD,mC,Len), scale))
  X1 <- detagene2_scaled[, c("SSD","LSD","mC","Len","Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(SSD=as.numeric(X1$SSD), LSD=as.numeric(X1$LSD),
                         mC=as.numeric(X1$mC), Len=as.numeric(X1$Len), Sub=factor(X1$Sub))
  train_idx <- createDataPartition(y_clean1, p=0.8, list=FALSE)
  data_train <- data.frame(X_clean1[train_idx,], Gene_mutation = y_clean1[train_idx])
  data_train <- data_train %>% distinct(SSD,LSD,mC,Len,Sub,Gene_mutation,.keep_all=TRUE)
  
  set.seed(123)
  rf_model <- ranger(Gene_mutation ~ ., data=data_train,
                     importance="permutation", num.trees=500,
                     mtry=3, min.node.size=10, oob.error=TRUE,
                     keep.inbag=TRUE, write.forest=TRUE)
  
  set.seed(123)
  sample_idx <- sample(seq_len(nrow(data_train)), size=min(5000, nrow(data_train)))
  X_sample <- data_train[sample_idx, c("SSD","LSD","mC","Len","Sub")]
  
  shap_values <- fastshap::explain(
    object = rf_model,
    X = X_sample,
    pred_wrapper = function(object,newdata) predict(object, data=newdata)$predictions,
    nsim = 500, adjust = TRUE
  )
  colnames(shap_values) <- gsub("^SHAP_SHAP_", "", colnames(shap_values))
  
  shap_abs_mean <- colMeans(abs(shap_values))
  shap_df <- data.frame(Feature=names(shap_abs_mean), Importance=shap_abs_mean) %>% arrange(desc(Importance))
  write.csv(shap_df, file.path(out_dir,"shap_importance.csv"), row.names=FALSE)
  
  p_bar <- ggplot(shap_df, aes(x=reorder(Feature, Importance), y=Importance, fill=Importance)) +
    geom_col(width=0.7) + coord_flip() +
    scale_fill_gradient(low="blue", high="red", guide="none") +
    labs(x=NULL, y="Mean |SHAP|") + theme_minimal(base_size=14)
  
  shap_long <- shap_values %>% as.data.frame() %>% mutate(.row=seq_len(n())) %>%
    tidyr::pivot_longer(-.row, names_to="Feature", values_to="SHAP") %>%
    left_join(
      X_sample %>% mutate(.row=seq_len(n())) %>% tidyr::pivot_longer(-.row, names_to="Feature", values_to="Value"),
      by=c(".row","Feature")
    )
  
  p_beeswarm <- ggplot(shap_long, aes(x=reorder(Feature, SHAP, FUN=function(x) median(abs(x))), y=SHAP, color=Value)) +
    geom_quasirandom(alpha=0.6, width=0.25, size=0.9, groupOnX=TRUE) +
    scale_color_gradient(low="blue", high="red", name="Feature Value") +
    coord_flip() + labs(x="Features", y="SHAP value", title="SHAP summary") + theme_minimal(base_size=14)
  
  p_final <- p_beeswarm | p_bar
  ggsave(file.path(out_dir,"shap_summary.png"), p_final, width=12, height=8, dpi=600, bg="white")
}

##3. rf_interaction_pipeline (交互作用分析)
run_rf_interaction <- function(data_file, out_dir) {
  dir.create(out_dir, showWarnings = FALSE)
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD,LSD,mC,Len), scale))
  X1 <- detagene2_scaled[, c("SSD","LSD","mC","Len","Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(SSD=as.numeric(X1$SSD), LSD=as.numeric(X1$LSD), mC=as.numeric(X1$mC), Len=as.numeric(X1$Len), Sub=factor(X1$Sub))
  
  train_idx <- createDataPartition(y_clean1, p=0.8, list=FALSE)
  data_train <- data.frame(X_clean1[train_idx,], Gene_mutation = y_clean1[train_idx])
  data_train <- data_train %>% distinct(SSD,LSD,mC,Len,Sub,Gene_mutation,.keep_all=TRUE)
  
  set.seed(123)
  rf_model <- ranger(Gene_mutation ~ ., data=data_train, importance="permutation", num.trees=500, mtry=3, min.node.size=10, oob.error=TRUE, keep.inbag=TRUE, write.forest=TRUE)
  
  predictor <- Predictor$new(model=rf_model, data=data_train[,c("SSD","LSD","mC","Len","Sub")], y=data_train$Gene_mutation)
  interaction_mc <- Interaction$new(predictor, feature="mC")
  write.csv(interaction_mc$results, file.path(out_dir,"interaction_results.csv"), row.names=FALSE)
}

## 4. permutation_test (置换检验) 
run_permutation_test <- function(data_file, n_perm, grid_points, num_cores, output_file) {
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD, LSD, mC, Len), scale))
  X1 <- detagene2_scaled[, c("SSD", "LSD", "mC", "Len", "Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(SSD=as.numeric(X1$SSD), LSD=as.numeric(X1$LSD), mC=as.numeric(X1$mC), Len=as.numeric(X1$Len), Sub=factor(X1$Sub))
  
  train_idx <- createDataPartition(y_clean1, p = 0.8, list = FALSE)
  data_train1 <- data.frame(X_clean1[train_idx, ], Gene_mutation = y_clean1[train_idx])
  data_train1 <- data_train1 %>% distinct(SSD,LSD,mC,Len,Sub,Gene_mutation,.keep_all=TRUE)
  
  set.seed(123)
  rf_model <- ranger(Gene_mutation ~ SSD+LSD+mC+Len+Sub, data=data_train1, importance="permutation", num.trees=500, mtry=3, min.node.size=10, oob.error=TRUE, keep.inbag=TRUE, write.forest=TRUE, num.threads=num_cores)
  
  registerDoParallel(num_cores)
  perm_h <- foreach(i = 1:n_perm, .combine = c, .packages = c("ranger","pdp")) %dopar% {
    data_perm <- data_train1
    data_perm$LSD <- sample(data_perm$LSD)
    rf_perm <- tryCatch({
      ranger(Gene_mutation ~ ., data=data_perm, num.trees=500, mtry=3, min.node.size=10)
    }, error=function(e) NULL)
    if (is.null(rf_perm)) return(NA)
    pd_perm <- pdp::partial(rf_perm, pred.var=c("LSD","mC"), grid.resolution=grid_points)
    if (nrow(pd_perm)<5) return(NA)
    var_total <- var(pd_perm$yhat, na.rm=TRUE)
    var_LSD <- var(pd_perm$yhat[pd_perm$mC == median(pd_perm$mC)], na.rm=TRUE)
    var_mc <- var(pd_perm$yhat[pd_perm$LSD == median(pd_perm$LSD)], na.rm=TRUE)
    if (var_total < 1e-8) 0 else (var_total - (var_LSD + var_mc)) / var_total
  }
  write.csv(perm_h, output_file)
  cat("置换检验结果已保存至:", output_file, "\n")
}

##主入口 
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  cat("用法: Rscript rf_master_pipeline.R <mode> <data_file.csv> <output_dir> [extra]\n")
  quit()
}
mode <- args[1]; data_file <- args[2]; out_dir <- args[3]

if (mode=="rf_pipeline") {
  run_rf_pipeline(data_file, out_dir)
} else if (mode=="rf_shap") {
  run_rf_shap(data_file, out_dir)
} else if (mode=="rf_interaction") {
  run_rf_interaction(data_file, out_dir)
} else if (mode=="permutation_test") {
  if (length(args)<6) { stop("permutation_test 模式需要: <n_perm> <grid_points> <num_cores> <output_file>") }
  run_permutation_test(data_file, as.integer(args[4]), as.integer(args[5]), as.integer(args[6]), args[7])
} else {
  stop("未知模式:", mode)
}
