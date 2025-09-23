# === 1. 加载必要的包 ===
library(dplyr)
library(pdp)
library(fastshap)
library(ranger)
library(caret)
library(doParallel)
library(foreach)

# === 2. 定义主要函数: run_permutation_test ===
run_permutation_test <- function(data_file, n_perm = 1000, grid_points = NULL, num_cores = 7, output_file = "perm_results.csv") {

  # 读取数据
  detagene2 <- read.csv(data_file)

  # 标准化数据
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD, LSD, mC, Len), scale))
  X1 <- detagene2_scaled[, c("SSD", "LSD", "mC", "Len", "Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  
  # Winsorize 极端值
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(
    SSD = as.numeric(X1$SSD),
    LSD = as.numeric(X1$LSD),
    mC = as.numeric(X1$mC),
    Len = as.numeric(X1$Len),
    Sub = factor(X1$Sub)
  )
  
  # 分层抽样平衡数据
  train_idx <- createDataPartition(y_clean1, p = 0.8, list = FALSE)
  X_train1 <- X_clean1[train_idx, ]
  y_train1 <- y_clean1[train_idx]
  epsilon <- 1e-16
  y_shift1 <- y_train1 + epsilon  # 避免 0
  data_train1 <- data.frame(X_train1, Gene_mutation = y_shift1)
  
  # 测试集
  X_ver1 <- X_clean1[-train_idx, ]
  y_ver1 <- y_clean1[-train_idx]
  y_shift_ver1 <- y_ver1 + epsilon
  data_ver1 <- data.frame(X_ver1, Gene_mutation = y_shift_ver1)
  
  # 删除重复样本（若为技术重复）
  dup_idx <- duplicated(data_train1) | duplicated(data_train1, fromLast = TRUE)
  cat("重复样本数量:", sum(dup_idx), "/", nrow(data_train1))
  data_train1 <- data_train1 %>% distinct(SSD, LSD, mC, Len, Sub, Gene_mutation, .keep_all = TRUE)

  # === 3. 随机森林建模 ===
  set.seed(123)
  rf_model <- ranger(
    formula = Gene_mutation ~ SSD + LSD + mC + Len + Sub,
    data = data_train1,
    importance = "permutation",
    num.trees = 500,
    mtry = 3,
    min.node.size = 10,
    oob.error = TRUE,
    keep.inbag = TRUE,
    write.forest = TRUE,
    num.threads = num_cores,
    verbose = FALSE
  )
  
  # === 4. 定义 H-statistic 计算 ===
  h_statistic <- function(pd_data) {
    var_total <- var(pd_data$yhat, na.rm = TRUE)
    var_LSD <- var(pd_data$yhat[pd_data$mC == median(pd_data$mC)], na.rm = TRUE)
    var_mc <- var(pd_data$yhat[pd_data$LSD == median(pd_data$LSD)], na.rm = TRUE)
    
    if (var_total < 1e-8) {  # 总方差接近零时返回0
      return(0)
    } else {
      (var_total - (var_LSD + var_mc)) / var_total
    }
  }

  # 初始化网格点
  if (is.null(grid_points)) {
    grid_points <- expand.grid(
      LSD = seq(min(data_train1$LSD), max(data_train1$LSD), length.out = 10),
      mC = seq(min(data_train1$mC), max(data_train1$mC), length.out = 10)
    )
  }

  # === 5. 进行置换检验 ===
  registerDoParallel(num_cores)
  
  perm_h <- foreach(i = 1:n_perm, .combine = c, .packages = c("ranger", "pdp")) %dopar% {
    data_perm <- data_train1
    data_perm$LSD <- sample(data_perm$LSD)
    
    rf_perm <- tryCatch({
      do.call(ranger, list(formula = Gene_mutation ~ ., data = data_perm, num.trees = 500, mtry = 3, min.node.size = 10, seed = 123))
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(rf_perm)) {
      return(NA)
    }
    
    pd_perm <- pdp::partial(rf_perm, pred.var = c("LSD", "mC"), pred.grid = grid_points)
    pd_perm <- pd_perm[complete.cases(pd_perm$yhat), ]
    
    if (nrow(pd_perm) < 5) {
      return(NA)
    }
    
    h_statistic(pd_perm)
  }

  # === 6. 保存结果 ===
  write.csv(perm_h, output_file)
  cat("置换检验结果已保存至:", output_file, "\n")
}

# === 7. 执行主脚本 ===
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]  # 数据文件路径
output_file <- args[2]  # 结果保存路径
n_perm <- as.integer(args[3])  # 置换次数
grid_points <- as.integer(args[4])  # 网格点
num_cores <- as.integer(args[5])  # 并行核心数

# 运行置换检验
run_permutation_test(data_file, n_perm, grid_points, num_cores, output_file)
