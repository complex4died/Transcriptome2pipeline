my_theme <- theme_bw(base_size = 14) +
  theme(
    # 图例
    legend.position  = "right",
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 14),
    # 坐标轴线与刻度
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    axis.line        = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    # 坐标轴文字
    axis.text.x      = element_text(angle = 0, size = 12),
    axis.text.y      = element_text(size = 12),
    axis.title.x     = element_text(size = 14),
    axis.title.y     = element_text(size = 14),
    # 标题
    plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 12, hjust = 0.5, color = "gray40"),
    # 分面标签
    strip.text       = element_text(size = 12),
    # 网格
    panel.grid       = element_blank(),
    # 全局字体
    text             = element_text(family = "sans")
  )




library(ggplot2)
library(ggrepel)
library(dplyr)

# ==============================================================================
# 火山图绘制函数 plot_volcano()
# ==============================================================================
#
#' @param deg_data   差异分析结果数据框（必须含 log2FC / pvalue / padj / Change 列）
#' @param pval_type  Y 轴使用的 P 值类型："padj"（默认）或 "pvalue"
#' @param pval_cutoff  显著性虚线阈值，默认 0.05
#' @param logfc_cutoff log2FC 虚线阈值，默认 1.5
#' @param colors     三元素命名向量，分别对应 Up / Stable / Down，可自定义
#' @param symmetric  X 轴是否对称（强制以最大绝对值为边界），默认 TRUE
#' @param top_n      是否自动标注 Top N 个显著基因，0 为不标注，默认 10
#' @param title      图标题，默认 "Volcano Plot"
#' @return ggplot 对象

plot_volcano <- function(deg_data,
                         pval_type    = c("padj", "pvalue"),
                         pval_cutoff  = 0.05,
                         logfc_cutoff = 1.5,
                         colors       = c("Up"     = "#ED0000B2",
                                          "Stable" = "grey60",
                                          "Down"   = "#00468BB2"),
                         symmetric    = TRUE,
                         top_n        = 10,
                         title        = "Volcano Plot") {
  
  pval_type <- match.arg(pval_type)
  
  # 1. 确认必要列存在
  required_cols <- c("log2FC", "pvalue", "padj", "Change")
  missing_cols  <- setdiff(required_cols, colnames(deg_data))
  if (length(missing_cols) > 0) {
    stop("数据框缺少必要列：", paste(missing_cols, collapse = ", "))
  }
  
  # 2. 动态指定 Y 轴数据列和标签
  y_col   <- pval_type
  y_label <- if (pval_type == "padj") "-log10(P.adjust)" else "-log10(P.value)"
  
  # 3. 过滤掉 Y 轴 NA（防止作图时报错）
  plot_data <- deg_data %>%
    dplyr::filter(!is.na(.data[[y_col]]) & !is.na(log2FC))
  
  # 4. X 轴边界：是否对称
  x_abs_max <- max(abs(plot_data$log2FC), na.rm = TRUE) * 1.05
  x_limits  <- if (symmetric) c(-x_abs_max, x_abs_max) else NULL
  
  # 5. Top 基因标注（Up 和 Down 分别取 top_n/2）
  label_data <- NULL
  if (top_n > 0 && "symbol" %in% colnames(plot_data)) {
    n_each <- ceiling(top_n / 2)
    
    top_up <- plot_data %>%
      dplyr::filter(Change == "Up") %>%
      dplyr::arrange(.data[[y_col]]) %>%   # padj/pvalue 越小越显著
      dplyr::slice_head(n = n_each)
    
    top_down <- plot_data %>%
      dplyr::filter(Change == "Down") %>%
      dplyr::arrange(.data[[y_col]]) %>%
      dplyr::slice_head(n = n_each)
    
    label_data <- dplyr::bind_rows(top_up, top_down)
  }
  
  # 6. 构建 ggplot 核心
  p <- ggplot(data = plot_data,
              aes(x = log2FC,
                  y = -log10(.data[[y_col]]),
                  fill = Change)) +
    
    geom_point(size   = 3,
               alpha  = 0.5,
               stroke = 0.2,
               shape  = 21,
               colour = "white") +
    
    scale_fill_manual(values = colors,
                      breaks = c("Up", "Stable", "Down")) +   # 固定图例顺序
    
    # X 轴：格式固定（log2FC），Y 轴字体格式固定（斜体 log10）
    labs(x     = expression(log[2] ~ "(Fold Change)"),
         y     = bquote(-log[10] ~ .(paste0("(", ifelse(pval_type == "padj", "P.adjust", "P.value"), ")"))),
         fill  = "Regulation",
         title = title) +
    
    # 阈值虚线
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
               linetype   = 2,
               color      = "gray30",
               linewidth  = 0.25) +
    geom_hline(yintercept = -log10(pval_cutoff),
               linetype   = 2,
               color      = "gray30",
               linewidth  = 0.25) +
    
    # 对称开关
    (if (symmetric) scale_x_continuous(limits = x_limits) else NULL) +
    
    my_theme
  
  # 7. 基因标注层（如有）
  if (!is.null(label_data) && nrow(label_data) > 0) {
    p <- p +
      geom_text_repel(
        data            = label_data %>% dplyr::filter(Change == "Up"),
        aes(label       = symbol),
        box.padding     = 0.5,
        max.overlaps    = Inf,
        size            = 3,
        fontface        = "italic",
        nudge_x         = 0.1,
        show.legend     = FALSE
      ) +
      geom_text_repel(
        data            = label_data %>% dplyr::filter(Change == "Down"),
        aes(label       = symbol),
        box.padding     = 0.5,
        max.overlaps    = Inf,
        size            = 3,
        fontface        = "italic",
        nudge_x         = -0.1,
        show.legend     = FALSE
      )
  }
  
  return(p)
}