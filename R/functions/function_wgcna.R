# ==============================================================================
# WGCNA 加权基因共表达网络分析统一模块
# 统一约定：
#   - 输入为 log 处理后的表达矩阵 (行名 = gene symbol, 列名 = 样本 ID)
#   - 返回 list，包含过滤后表达矩阵、软阈值、模块颜色、模块特征基因 (MEs) 及
#     模块-性状相关性矩阵；绘图函数 plot_module_trait 单独调用
# 说明：WGCNA 会 mask 部分基础函数，故本文件一律用 WGCNA:: 前缀，不 attach 包
# ==============================================================================

# 内部工具：将样本信息表转换为可用于相关性分析的数值性状矩阵
#  - 分组列 -> 0/1 (1 = case)
#  - 数值列 -> 原样
#  - 二分类字符/因子列 -> 0/1
#  - 多分类列 -> 哑变量 (one-hot)
build_traits <- function(meta, sample_id_col, group_col, case, trait_cols) {
  out <- data.frame(row.names = meta[[sample_id_col]])
  out[[case]] <- as.numeric(meta[[group_col]] == case)

  for (tc in setdiff(trait_cols, group_col)) {
    x <- meta[[tc]]
    if (is.numeric(x)) {
      out[[tc]] <- as.numeric(x)
    } else {
      x <- factor(x)
      lv <- levels(x)
      if (length(lv) == 2) {
        out[[tc]] <- as.numeric(x) - 1
      } else {
        for (l in lv) out[[paste0(tc, "_", l)]] <- as.numeric(x == l)
      }
    }
  }
  out
}

#' WGCNA 加权共表达网络构建 + 模块识别 + 模块-性状相关性
#'
#' @param exp_mat log 处理后的表达矩阵 (行名 = gene symbol, 列名 = 样本 ID)
#' @param sample_info 样本信息表 (需包含样本 ID 列、分组列及性状列)
#' @param group_col 分组列名, 默认 "Groups"
#' @param case 实验组名称, 默认 "sepsis"
#' @param control 对照组名称, 默认 "healthy"
#' @param sample_id_col 样本 ID 列名, 默认 "Sample_ID"
#' @param trait_cols 附加性状列 (连续或分类均可, 自动转换), 默认 c("Age", "Gender")
#' @param n_var_genes 保留方差最大的前 N 个基因, 默认 8000; NULL 表示使用全部基因
#' @param network_type 网络类型 "signed" / "unsigned" / "signed hybrid", 默认 "signed"
#' @param soft_power 软阈值; NULL 表示按无标度拓扑自动选择, 默认 NULL
#' @param power_vector 软阈值候选网格, 默认 1:20
#' @param n_threads WGCNA 并行线程数, 默认 4
#' @param rsquared_cut 无标度拓扑拟合 R^2 阈值, 默认 0.85
#' @param min_module_size 最小模块基因数, 默认 30
#' @param merge_cut_height 模块合并树高阈值 (越小越倾向多模块), 默认 0.25
#' @param max_block_size 单块最大基因数 (大于此值自动分块), 默认 20000
#' @param numeric_labels 是否使用数值模块标签 (TRUE) 或颜色标签 (FALSE), 默认 TRUE
#' @param save_toms 是否保存 TOM 矩阵到磁盘 (后续做 hub gene / Cytoscape 需要), 默认 FALSE
#' @param tom_file_base TOM 保存文件名前缀, 默认 "wgcnaTOM"
#' @param verbose 详细输出等级 (0 = 静默, 3 = 详细), 默认 0
#'
#' @return list, 包含 dat_expr / soft_power / sft / net / module_colors /
#'         module_labels / MEs / traits / module_trait_cor / module_trait_pval
run_wgcna <- function(exp_mat,
                      sample_info,
                      group_col        = "Groups",
                      case             = "sepsis",
                      control          = "healthy",
                      sample_id_col    = "Sample_ID",
                      trait_cols       = c("Age", "Gender"),
                      n_var_genes      = 8000,
                      network_type     = "signed",
                      soft_power       = NULL,
                      power_vector     = 1:20,
                      n_threads        = 4,
                      rsquared_cut     = 0.85,
                      min_module_size  = 30,
                      merge_cut_height = 0.25,
                      max_block_size   = 20000,
                      numeric_labels   = TRUE,
                      save_toms        = FALSE,
                      tom_file_base    = "wgcnaTOM",
                      verbose          = 0) {

  if (!requireNamespace("WGCNA", quietly = TRUE))
    stop("缺少 WGCNA 包，请先运行 BiocManager::install('WGCNA')", call. = FALSE)

  # 关键：必须 attach WGCNA（不能只用 WGCNA:: 前缀）。
  # WGCNA 内部通过 match.fun("cor") 查找 cor 函数；若不 attach，
  # S4Vectors/IRanges 会把 stats::cor 变成 S4 泛型并拦截该查找，
  # 导致报错 "unused arguments (weights.x, weights.y, cosine)"。
  suppressPackageStartupMessages(library(WGCNA))

  # 开启多线程（默认 4 核）
  WGCNA::allowWGCNAThreads(nThreads = n_threads)

  # 1. 转成 WGCNA 要求的 样本 x 基因 布局
  dat_expr <- as.data.frame(t(as.matrix(exp_mat)))

  # 2. 过滤低方差基因，保留 top n_var_genes
  if (!is.null(n_var_genes) && n_var_genes < ncol(dat_expr)) {
    gene_var <- apply(dat_expr, 2, stats::var, na.rm = TRUE)
    keep <- order(gene_var, decreasing = TRUE)[seq_len(n_var_genes)]
    dat_expr <- dat_expr[, keep, drop = FALSE]
  }

  # 3. 质控：剔除缺失过多的基因 / 样本
  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = verbose)
  dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

  # 4. 软阈值选择
  if (is.null(soft_power)) {
    sft <- WGCNA::pickSoftThreshold(
      dat_expr,
      powerVector = power_vector,
      networkType = network_type,
      RsquaredCut = rsquared_cut,
      verbose = verbose
    )
    soft_power <- sft$powerEstimate
    if (is.na(soft_power)) {
      soft_power <- 6
      message(sprintf(">>> 未找到满足 R^2 >= %.2f 的软阈值，回退使用 soft power = %d", rsquared_cut, soft_power))
    } else if (soft_power <= 2) {
      message(sprintf(">>> 自动选出的软阈值 soft power = %d 偏低，网络可能过密；建议结合 mean connectivity 手动指定更大的 soft_power（如 6-8）。", soft_power))
    }
  } else {
    sft <- NULL
  }

  # 5. 网络构建 + 模块识别
  net <- WGCNA::blockwiseModules(
    dat_expr,
    power              = soft_power,
    networkType        = network_type,
    TOMType            = network_type,
    minModuleSize      = min_module_size,
    reassignThreshold  = 0,
    mergeCutHeight     = merge_cut_height,
    numericLabels      = numeric_labels,
    pamRespectsDendro  = FALSE,
    maxBlockSize       = max_block_size,
    saveTOMs           = save_toms,
    saveTOMFileBase    = tom_file_base,
    verbose            = verbose
  )

  module_colors <- WGCNA::labels2colors(net$colors)
  MEs           <- WGCNA::orderMEs(net$MEs)

  # 6. 构建性状矩阵并做模块-性状相关性
  traits <- build_traits(sample_info, sample_id_col, group_col, case, trait_cols)
  traits <- traits[rownames(dat_expr), , drop = FALSE]

  module_trait_cor  <- stats::cor(MEs, traits, use = "pairwise.complete.obs")
  module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nSamples = nrow(dat_expr))

  list(
    dat_expr          = dat_expr,
    soft_power        = soft_power,
    sft               = sft,
    net               = net,
    module_colors     = module_colors,
    module_labels     = net$colors,
    MEs               = MEs,
    traits            = traits,
    module_trait_cor  = module_trait_cor,
    module_trait_pval = module_trait_pval
  )
}

#' 模块-性状相关性热图 (基于 WGCNA::labeledHeatmap)
#'
#' @param wgcna_res run_wgcna 返回的 list
#' @param zlim 颜色映射范围, 默认 c(-1, 1)
#' @param cex_text 热图内文字大小, 默认 0.7
#' @param title 图标题, 默认 "Module-trait relationships"
#'
#' @return 直接绘制到当前图形设备, 不可见返回相关性矩阵
plot_module_trait <- function(wgcna_res,
                              zlim     = c(-1, 1),
                              cex_text = 0.7,
                              title    = "Module-trait relationships") {

  cor_mat <- wgcna_res$module_trait_cor
  pval    <- wgcna_res$module_trait_pval

  WGCNA::labeledHeatmap(
    Matrix       = cor_mat,
    xLabels      = colnames(cor_mat),
    yLabels      = rownames(cor_mat),
    ySymbols     = rownames(cor_mat),
    colorLabels  = FALSE,
    colors       = WGCNA::blueWhiteRed(50),
    textMatrix   = paste0(signif(cor_mat, 2), "\n(", signif(pval, 1), ")"),
    setStdMargins = FALSE,
    cex.text     = cex_text,
    zlim         = zlim,
    main         = title
  )

  invisible(cor_mat)
}
