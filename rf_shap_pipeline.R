#!/usr/bin/env Rscript
# ===============================================
# 机器学习shap重要性分析脚本
# 作者: 杜伟轩
# 功能:
#   1. 支持随机森林模型
#   2. 支持shap解释模型
#   3. 支持自动生成可视化
# ===============================================
##使用Rscript rf_shap_pipeline.R ALLdistanceHC9.csv output_dir/
suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(ranger)
  library(Boruta)
  library(fastshap)
  library(shapviz)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(ggbeeswarm)
})

#主函数 
run_rf_shap <- function(data_file, out_dir = "results_out") {
  dir.create(out_dir, showWarnings = FALSE)

  # 数据读取与处理
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD,LSD,mC,Len), scale))
  X1 <- detagene2_scaled[, c("SSD","LSD","mC","Len","Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(SSD=as.numeric(X1$SSD),
                         LSD=as.numeric(X1$LSD),
                         mC=as.numeric(X1$mC),
                         Len=as.numeric(X1$Len),
                         Sub=factor(X1$Sub))

  # 分层抽样
  train_idx <- createDataPartition(y_clean1, p=0.8, list=FALSE)
  data_train <- data.frame(X_clean1[train_idx,], Gene_mutation = y_clean1[train_idx])
  data_train <- data_train %>% distinct(SSD,LSD,mC,Len,Sub,Gene_mutation,.keep_all=TRUE)

  # 模型训练 
  set.seed(123)
  rf_model <- ranger(Gene_mutation ~ ., data=data_train,
                     importance="permutation", num.trees=500,
                     mtry=3, min.node.size=10, oob.error=TRUE,
                     keep.inbag=TRUE, write.forest=TRUE)

  # 计算 SHAP 值 
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
  shap_df <- data.frame(Feature=names(shap_abs_mean),
                        Importance=shap_abs_mean) %>% arrange(desc(Importance))
  write.csv(shap_df, file.path(out_dir,"shap_importance.csv"), row.names=FALSE)

  # 绘图 
  # 条形图（红蓝渐变）
  p_bar <- ggplot(shap_df, aes(x=reorder(Feature, Importance), y=Importance, fill=Importance)) +
    geom_col(width=0.7) + coord_flip() +
    scale_fill_gradient(low="blue", high="red", guide="none") +
    labs(x=NULL, y="Mean |SHAP|") + theme_minimal(base_size=14)

  # Beeswarm 图
  shap_long <- shap_values %>% as.data.frame() %>%
    mutate(.row=seq_len(n())) %>%
    tidyr::pivot_longer(-.row, names_to="Feature", values_to="SHAP") %>%
    left_join(
      X_sample %>% mutate(.row=seq_len(n())) %>%
        tidyr::pivot_longer(-.row, names_to="Feature", values_to="Value"),
      by=c(".row","Feature")
    )

  p_beeswarm <- ggplot(shap_long, aes(x=reorder(Feature, SHAP, FUN=function(x) median(abs(x))),
                                      y=SHAP, color=Value)) +
    geom_quasirandom(alpha=0.6, width=0.25, size=0.9, groupOnX=TRUE) +
    scale_color_gradient(low="blue", high="red", name="Feature Value") +
    coord_flip() + labs(x="Features", y="SHAP value", title="SHAP summary") +
    theme_minimal(base_size=14)

  # 拼图并保存
  p_final <- p_beeswarm | p_bar
  ggsave(file.path(out_dir,"shap_summary.png"), p_final, width=12, height=8, dpi=600, bg="white")
}

# 脚本执行入口 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript rf_shap_pipeline.R <data_file.csv> <output_dir>\n")
} else {
  run_rf_shap(args[1], args[2])
}
