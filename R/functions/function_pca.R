library(ggplot2)
library(dplyr)

#' 绘制 PCA 主成分分析散点图 (沿用 Volcano Plot 视觉风格)
#'
#' @param exp_mat log 处理后的表达矩阵 (行名为基因, 列名为样本 ID, 如 exp_cpm_log1p)
#' @param sample_info 样本信息表 (需包含 Sample_ID 和分组列)
#' @param group_col 分组列名, 默认 "Groups"
#' @param case 实验组名称, 默认 "sepsis"
#' @param control 对照组名称, 默认 "healthy"
#' @param group_colors 组别颜色映射命名向量, 默认按 control/case 映射; 也可直接传入具体组名如 c("healthy" = "#00468BB2", "sepsis" = "#ED0000B2")
#' @param top_n_var 可选, 仅使用方差最大的前 N 个基因计算 PCA; 默认 NULL 表示全部有效基因
#' @param scale_genes 是否对基因做标准化 (Z-score), 默认 TRUE
#' @param show_ellipse 是否绘制 95% 置信椭圆, 默认 TRUE
#' @param title 图表标题, 默认 "PCA Plot"
#'
#' @return ggplot 对象
plot_pca <- function(exp_mat,
                     sample_info,
                     group_col    = "Groups",
                     case         = "sepsis",
                     control      = "healthy",
                     group_colors = c("control" = "#00468BB2",
                                      "case"    = "#ED0000B2"),
                     top_n_var    = NULL,
                     scale_genes  = TRUE,
                     show_ellipse = TRUE,
                     title        = "PCA Plot") {

  # 1. 对齐样本
  common_samples <- intersect(colnames(exp_mat), sample_info$Sample_ID)
  if (length(common_samples) < 3) stop("有效匹配样本数不足 3 个，无法计算 PCA！")

  mat <- as.matrix(exp_mat[, common_samples, drop = FALSE])
  meta_subset <- sample_info[match(common_samples, sample_info$Sample_ID), ]

  # 2. 过滤零方差 / NA 基因，避免 prcomp(scale. = TRUE) 报错
  gene_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  valid_idx <- !is.na(gene_vars) & gene_vars > 0
  mat <- mat[valid_idx, , drop = FALSE]
  gene_vars <- gene_vars[valid_idx]

  if (!is.null(top_n_var) && top_n_var < nrow(mat)) {
    top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(top_n_var)]
    mat <- mat[top_idx, , drop = FALSE]
  }

  # 3. 计算 PCA
  pca_res <- stats::prcomp(t(mat), center = TRUE, scale. = scale_genes)
  var_pct <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

  pca_df <- data.frame(
    Sample_ID = common_samples,
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Group = factor(meta_subset[[group_col]], levels = c(control, case)),
    stringsAsFactors = FALSE
  )

  # 4. 组别颜色映射 (支持默认 "control"/"case" 占位键名、具体组名命名向量或无名向量)
  if (is.null(names(group_colors))) {
    group_colors <- stats::setNames(group_colors, c(control, case))
  } else if (all(c("control", "case") %in% names(group_colors))) {
    group_colors <- stats::setNames(
      unname(group_colors[c("control", "case")]),
      c(control, case)
    )
  }

  # 5. 绘图 (沿用 Volcano Plot 风格)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Group)) +
    geom_vline(xintercept = 0, linetype = 2, color = "gray30", linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray30", linewidth = 0.25) +
    geom_point(
      size = 3,
      alpha = 0.7,
      stroke = 0.2,
      shape = 21,
      colour = "white"
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(
      x = paste0("PC1 (", var_pct[1], "%)"),
      y = paste0("PC2 (", var_pct[2], "%)"),
      title = title
    )

  if (show_ellipse) {
    p <- p + stat_ellipse(
      aes(color = Group),
      level = 0.95,
      linetype = 2,
      linewidth = 0.5,
      show.legend = FALSE
    )
  }

  if (exists("my_theme", inherits = TRUE)) {
    p <- p + my_theme
  } else {
    p <- p + theme_classic(base_size = 14)
  }

  return(p)
}
