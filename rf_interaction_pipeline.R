#!/usr/bin/env Rscript
# ===============================================
# 机器学习单因素与多因素互作分析脚本
# 作者: 杜伟轩
# 功能:
#   1. 支持随机森林模型
#   2. 支持iml包解释模型
#   3. 支持自动生成可视化
# ===============================================
##使用Rscript rf_interaction_pipeline.R ALLdistanceHC9.csv output_dir/
suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(ranger)
  library(Boruta)
  library(car)
  library(pdp)
  library(iml)
  library(ggplot2)
  library(patchwork)
})

# 主函数
run_pipeline <- function(data_file, out_dir = "results_out") {
  dir.create(out_dir, showWarnings = FALSE)

  #数据读取
  detagene2 <- read.csv(data_file)
  detagene2_scaled <- detagene2 %>% mutate(across(c(SSD,LSD,mC,Len), scale))
  X1 <- detagene2_scaled[, c("SSD","LSD","mC","Len","Sub")]
  y1 <- detagene2_scaled$Gene_mutation
  y_clean1 <- ifelse(y1 > quantile(y1, 0.99), quantile(y1, 0.99), y1)
  X_clean1 <- data.frame(SSD=as.numeric(X1$SSD), LSD=as.numeric(X1$LSD),
                         mC=as.numeric(X1$mC), Len=as.numeric(X1$Len),
                         Sub=factor(X1$Sub))

  # 数据集划分
  train_idx <- createDataPartition(y_clean1, p=0.8, list=FALSE)
  data_train <- data.frame(X_clean1[train_idx,], Gene_mutation = y_clean1[train_idx])
  data_train <- data_train %>% distinct(SSD,LSD,mC,Len,Sub,Gene_mutation,.keep_all=TRUE)

  #  随机森林建模 
  set.seed(123)
  rf_model <- ranger(Gene_mutation ~ ., data=data_train,
                     importance="permutation", num.trees=500,
                     mtry=3, min.node.size=10, oob.error=TRUE,
                     keep.inbag=TRUE, write.forest=TRUE)

  # PDP 单变量
  pd_list <- list(
    LSD = partial(rf_model, pred.var="LSD", grid.resolution=30),
    mC  = partial(rf_model, pred.var="mC", grid.resolution=30),
    SSD = partial(rf_model, pred.var="SSD", grid.resolution=30),
    Len = partial(rf_model, pred.var="Len", grid.resolution=30)
  )

  # 绘制单变量效应
  plot_pdp <- function(pd, var, title) {
    ggplot(pd, aes_string(x=var, y="yhat")) +
      geom_point(alpha=0.6, size=2, aes(color=yhat)) +
      geom_smooth(method="loess", color="purple", fill="pink", se=TRUE) +
      scale_color_gradient(low="blue", high="red") +
      labs(title=title, y="Prediction", x=var) +
      theme_minimal(base_size=14)
  }

  p1 <- plot_pdp(pd_list$mC, "mC", "mC effect")
  p2 <- plot_pdp(pd_list$LSD, "LSD", "LSD effect")
  p3 <- plot_pdp(pd_list$SSD, "SSD", "SSD effect")
  p4 <- plot_pdp(pd_list$Len, "Len", "Length effect")

  ggsave(file.path(out_dir,"pdp_features.png"), (p1/p2) | (p3/p4), width=12, height=8, dpi=600)

  # --- 交互作用分析 (iml) ---
  predictor <- Predictor$new(model=rf_model, data=data_train[,c("SSD","LSD","mC","Len","Sub")],
                             y=data_train$Gene_mutation)
  interaction_mc <- Interaction$new(predictor, feature="mC")
  write.csv(interaction_mc$results, file.path(out_dir,"interaction_results.csv"), row.names=FALSE)

  #  PDP 双变量交互
  pd_int <- list(
    LSD_mC = partial(rf_model, pred.var=c("LSD","mC"), grid.resolution=30),
    SSD_mC = partial(rf_model, pred.var=c("SSD","mC"), grid.resolution=30),
    Len_mC = partial(rf_model, pred.var=c("Len","mC"), grid.resolution=30)
  )

  plot_int <- function(pd, x, y, title, hval) {
    ggplot(pd, aes_string(x=x, y=y, fill="yhat")) +
      geom_tile(color="white", linewidth=0.3) +
      scale_fill_gradient(low="blue", high="red") +
      labs(title=title, subtitle=paste0("H = ", round(hval,3))) +
      theme_minimal(base_size=14) +
      coord_fixed(ratio=1)
  }

  hvals <- interaction_mc$results
  p_int1 <- plot_int(pd_int$LSD_mC,"LSD","mC","mC-LSD Interaction",
                     hvals$.interaction[hvals$.feature=="LSD:mC"])
  p_int2 <- plot_int(pd_int$SSD_mC,"SSD","mC","mC-SSD Interaction",
                     hvals$.interaction[hvals$.feature=="SSD:mC"])
  p_int3 <- plot_int(pd_int$Len_mC,"Len","mC","mC-Len Interaction",
                     hvals$.interaction[hvals$.feature=="Len:mC"])

  ggsave(file.path(out_dir,"interaction_plots.png"), p_int1 | p_int2 | p_int3, width=14, height=6, dpi=600)
}

# 脚本执行入口
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  cat("用法: Rscript rf_interaction_pipeline.R <data_file.csv> <output_dir>\n")
} else {
  run_pipeline(args[1], args[2])
}
