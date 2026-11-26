library(ComplexHeatmap)
library(circlize)
library(dplyr)

#' 绘制差异基因热图 (ComplexHeatmap 封装版)
#'
#' @param exp_mat 表达矩阵 (行名为基因 symbol, 列名为样本 ID)
#' @param deg_table 差异分析结果表 (需包含 symbol, Change, log2FC 列)
#' @param sample_info 样本信息表 (需包含 Sample_ID 和分组列)
#' @param group_col 分组列名, 默认 "Groups"
#' @param case 实验组名称
#' @param control 对照组名称
#' @param top_n 每组提取的 top 差异基因数量, 默认 50 (总共显示 2 * top_n 个)
#' @param group_colors 组别颜色映射命名向量, 例如 c(healthy = "#2166AC", sepsis = "#B2182B")
#' @param split_by_group 是否按分组切分列 (Control 在左, Case 在右, 保留组内聚类), 默认 TRUE
#' @param title 图表标题
#'
#' @return ComplexHeatmap 对象
plot_heatmap <- function(exp_mat, 
                         deg_table, 
                         sample_info, 
                         group_col = "Groups", 
                         case = "sepsis", 
                         control = "healthy", 
                         top_n = 50, 
                         group_colors = NULL, 
                         split_by_group = TRUE,
                         title = "Differential Expression Heatmap") {
  
  # 1. 筛选 Top 基因
  # 分别提取 Up 和 Down 基因
  top_up   <- deg_table |> 
    dplyr::filter(Change == "Up") |> 
    dplyr::arrange(desc(abs(log2FC))) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::pull(symbol)
  
  top_down <- deg_table |> 
    dplyr::filter(Change == "Down") |> 
    dplyr::arrange(desc(abs(log2FC))) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::pull(symbol)
  
  plot_genes <- c(top_up, top_down)
  
  # 2. 提取并标准化数据
  # 确保基因在矩阵中存在
  valid_genes <- intersect(plot_genes, rownames(exp_mat))
  if (length(valid_genes) == 0) stop("没有筛选出的基因在表达矩阵中找到！")
  
  # 提取数据
  plot_mat <- exp_mat[valid_genes, ]
  
  # Z-score 标准化 (Scale by row)
  # t(scale(t(mat))) 是 R 中对行做 scale 的标准写法
  plot_mat_scaled <- t(scale(t(as.matrix(plot_mat))))
  
  # 3. 颜色映射定义
  # 定义 -2 (蓝) -> 0 (白) -> 2 (红)
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
  
  # 4. 样本对齐与列注释
  common_samples <- intersect(colnames(plot_mat), sample_info$Sample_ID)
  plot_mat_scaled <- plot_mat_scaled[, common_samples]
  meta_subset <- sample_info[match(common_samples, sample_info$Sample_ID), ]
  
  # 设置分组颜色
  if (is.null(group_colors)) {
    group_colors <- stats::setNames(c("#2166AC", "#B2182B"), c(control, case))
  }
  
  anno_col <- ComplexHeatmap::HeatmapAnnotation(
    Group = meta_subset[[group_col]],
    col = list(Group = group_colors)
  )
  
  col_split <- if (split_by_group) {
    factor(meta_subset[[group_col]], levels = c(control, case))
  } else {
    NULL
  }
  
  # 5. 绘制热图
  ht <- ComplexHeatmap::Heatmap(
    plot_mat_scaled, 
    name = "Z-score", 
    col = col_fun, 
    
    # 行设置
    cluster_rows = TRUE, 
    show_row_names = (length(valid_genes) <= 100), 
    row_names_gp = gpar(fontsize = 8), 
    
    # 列设置
    cluster_columns = TRUE, 
    column_split = col_split,
    cluster_column_slices = FALSE,
    show_column_names = FALSE, 
    top_annotation = anno_col, 
    
    # 样式
    border = TRUE, 
    column_title = title,
    row_title = "Top DEGs",
    column_gap = unit(2, "mm")
  )

  
  return(ht)
}
