# ==============================================================================
# 富集分析统一模块：GO / KEGG × ORA / GSEA（基于 clusterProfiler）
# 统一约定：
#   - 所有函数输入的基因 / 背景均为 gene symbol
#   - ORA  输入显著基因列表 (symbol)；GSEA 输入全基因排序向量 (名字为 symbol)
#   - 均返回 list(table = 结果数据框, object = clusterProfiler 对象)
# 说明：clusterProfiler / AnnotationDbi 会 mask dplyr 的 select/rename/filter，
#       为避免污染全局命名空间，本文件一律用命名空间限定符调用，不 attach 包。
# ==============================================================================

# -----------------------------------------------------------------------------
# 0. 内部工具：物种 OrgDb 匹配 & symbol -> Entrez 映射（统一使用 clusterProfiler::bitr）
# -----------------------------------------------------------------------------
get_orgdb <- function(organism = "hsa") {
  pkg <- switch(
    tolower(organism),
    hsa = "org.Hs.eg.db",
    mmu = "org.Mm.eg.db",
    rno = "org.Rn.eg.db",
    stop("暂不支持的 organism: ", organism, "（当前支持 hsa / mmu / rno）", call. = FALSE)
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少物种注释包 ", pkg, "，请先运行 BiocManager::install('", pkg, "')", call. = FALSE)
  }
  getExportedValue(pkg, pkg)
}

symbol_to_entrez <- function(symbols, org_db = org.Hs.eg.db::org.Hs.eg.db) {
  symbols <- unique(trimws(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) stop("输入的基因列表为空！")

  map <- suppressMessages(suppressWarnings(
    clusterProfiler::bitr(symbols,
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org_db)
  ))
  entrez <- unique(map$ENTREZID[!is.na(map$ENTREZID)])
  if (length(entrez) == 0) stop("没有任何 symbol 成功映射到 Entrez ID，请检查物种或基因名格式！")
  message(sprintf(">>> [ID转换] %d 个 symbol 映射到 %d 个唯一 Entrez ID",
                  length(symbols), length(entrez)))
  entrez
}

# 将 symbol 命名的排序向量转换为 Entrez 命名（GSEA-KEGG 用）
rank_to_entrez <- function(ranked, org_db = org.Hs.eg.db::org.Hs.eg.db) {
  map <- suppressMessages(suppressWarnings(
    clusterProfiler::bitr(names(ranked),
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org_db)
  ))
  map <- map[!is.na(map$ENTREZID) & !is.na(map$SYMBOL), , drop = FALSE]
  if (nrow(map) == 0) stop("没有任何基因成功映射到 Entrez ID，无法运行 GSEA-KEGG！")
  map$score <- ranked[map$SYMBOL]

  # 多个 symbol 映射到同一 Entrez 时，保留 |score| 最大的一个
  map <- map[order(abs(map$score), decreasing = TRUE), , drop = FALSE]
  map <- map[!duplicated(map$ENTREZID), , drop = FALSE]

  out <- stats::setNames(map$score, map$ENTREZID)
  out <- out[order(out, decreasing = TRUE)]
  message(sprintf(">>> [ID转换] %d 个基因保留 %d 个唯一 Entrez ID 用于 GSEA",
                  length(ranked), length(out)))
  out
}

# 统一收尾 KEGG 结果：安全处理 NULL / 空结果，再转回可读 Gene Symbol
finalize_kegg_result <- function(res, org_db) {
  if (is.null(res)) {
    warning("KEGG 未检测到有效映射基因或无通路结果，返回空结果。")
    return(list(table = data.frame(), object = NULL))
  }
  df <- as.data.frame(res)
  if (nrow(df) > 0) {
    res <- clusterProfiler::setReadable(res, OrgDb = org_db, keyType = "ENTREZID")
    df  <- as.data.frame(res)
  }
  list(table = df, object = res)
}

# -----------------------------------------------------------------------------
# 1. 输入准备函数（从 DEG 结果 / 表达矩阵自动生成）
# -----------------------------------------------------------------------------

#' 提取显著差异基因列表（ORA 输入）
#' @param deg_table run_deg 返回的差异分析表（需含 symbol / Change 列）
#' @param direction "both" 合并 Up+Down；"up" 仅上调；"down" 仅下调
get_sig_genes <- function(deg_table, direction = c("both", "up", "down")) {
  direction <- match.arg(direction)
  keep <- switch(direction,
                 both = c("Up", "Down"),
                 up   = "Up",
                 down = "Down")
  deg_table$symbol[deg_table$Change %in% keep]
}

#' 生成全基因排序向量（GSEA 输入）
#' @param deg_table 差异分析表（需含 symbol / log2FC 列）
#' @return 命名数值向量，按 log2FC 降序，名字为 symbol
rank_genes <- function(deg_table) {
  ranked <- deg_table[!is.na(deg_table$log2FC), ]
  ranked <- ranked[!duplicated(ranked$symbol), ]
  out <- stats::setNames(ranked$log2FC, ranked$symbol)
  out[order(out, decreasing = TRUE)]
}

#' 提取背景基因集（ORA 的 universe）
#' @param exp_mat 表达矩阵（行名为 gene symbol）
get_universe <- function(exp_mat) {
  rownames(exp_mat)
}

# 内部：按 pval_type ("padj" / "pvalue") 统一过滤富集对象
filter_enrich_by_pval <- function(res, pval_type = c("padj", "pvalue"), pval_cutoff = 0.05) {
  pval_type <- match.arg(pval_type)
  if (is.null(res)) return(NULL)
  p_col <- if (pval_type == "padj") "p.adjust" else "pvalue"
  df <- res@result
  if (nrow(df) > 0 && p_col %in% colnames(df)) {
    res@result <- df[!is.na(df[[p_col]]) & df[[p_col]] <= pval_cutoff, , drop = FALSE]
  }
  if ("pvalueCutoff" %in% slotNames(res)) {
    res@pvalueCutoff <- pval_cutoff
  }
  res
}

# -----------------------------------------------------------------------------
# 2. ORA 富集（过表达分析）
# -----------------------------------------------------------------------------

#' GO 过表达分析（ORA）
#' @param genes 显著基因列表 (symbol)
#' @param universe 背景基因集 (symbol)，建议传入表达矩阵全部基因
#' @param ont GO 子本体："BP" / "MF" / "CC" / "ALL"
#' @param pval_type 筛选指标："padj"（默认，校正后 p 值）或 "pvalue"（原始 p 值）
enrich_go <- function(genes, universe = NULL, organism = "hsa", ont = "BP",
                      pval_cutoff = 0.05, qval_cutoff = 0.2,
                      pval_type = c("padj", "pvalue"),
                      p_adjust_method = "BH", min_gs_size = 10, max_gs_size = 500) {
  pval_type <- match.arg(pval_type)
  org_db    <- get_orgdb(organism)
  res <- clusterProfiler::enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org_db,
    keyType       = "SYMBOL",
    ont           = ont,
    pvalueCutoff  = if (pval_type == "pvalue") 1 else pval_cutoff,
    pAdjustMethod = p_adjust_method,
    qvalueCutoff  = if (pval_type == "pvalue") 1 else qval_cutoff,
    minGSSize     = min_gs_size,
    maxGSSize     = max_gs_size,
    readable      = TRUE
  )
  res <- filter_enrich_by_pval(res, pval_type = pval_type, pval_cutoff = pval_cutoff)
  if (is.null(res)) return(list(table = data.frame(), object = NULL))
  list(table = as.data.frame(res), object = res)
}

#' KEGG 过表达分析（ORA，基于 clusterProfiler::enrichKEGG）
#' @param genes 显著基因列表 (symbol)
#' @param universe 背景基因集 (symbol)
#' @param pval_type 筛选指标："padj"（默认）或 "pvalue"
enrich_kegg <- function(genes, universe = NULL, organism = "hsa",
                        pval_cutoff = 0.05, qval_cutoff = 0.2,
                        pval_type = c("padj", "pvalue"),
                        p_adjust_method = "BH", min_gs_size = 10, max_gs_size = 500,
                        use_internal_data = FALSE) {
  pval_type       <- match.arg(pval_type)
  org_db          <- get_orgdb(organism)
  gene_entrez     <- symbol_to_entrez(genes, org_db = org_db)
  universe_entrez <- if (!is.null(universe)) symbol_to_entrez(universe, org_db = org_db) else NULL

  if (is.null(getOption("clusterProfiler.download.method"))) {
    options(clusterProfiler.download.method = "curl")
  }

  res <- clusterProfiler::enrichKEGG(
    gene              = gene_entrez,
    universe          = universe_entrez,
    organism          = organism,
    keyType           = "ncbi-geneid",
    pvalueCutoff      = if (pval_type == "pvalue") 1 else pval_cutoff,
    pAdjustMethod     = p_adjust_method,
    qvalueCutoff      = if (pval_type == "pvalue") 1 else qval_cutoff,
    minGSSize         = min_gs_size,
    maxGSSize         = max_gs_size,
    use_internal_data = use_internal_data
  )
  res <- filter_enrich_by_pval(res, pval_type = pval_type, pval_cutoff = pval_cutoff)
  finalize_kegg_result(res, org_db = org_db)
}

# -----------------------------------------------------------------------------
# 3. GSEA 富集（基因集富集分析）
# -----------------------------------------------------------------------------

#' GO 基因集富集分析（GSEA）
#' @param ranked 全基因排序向量（名字为 symbol，按 log2FC 降序）
#' @param pval_type 筛选指标："padj"（默认）或 "pvalue"
gsea_go <- function(ranked, organism = "hsa", ont = "BP", pval_cutoff = 0.05,
                    pval_type = c("padj", "pvalue"),
                    p_adjust_method = "BH", min_gs_size = 10, max_gs_size = 500, ...) {
  pval_type <- match.arg(pval_type)
  org_db    <- get_orgdb(organism)
  res <- clusterProfiler::gseGO(
    geneList      = ranked,
    OrgDb         = org_db,
    keyType       = "SYMBOL",
    ont           = ont,
    pvalueCutoff  = if (pval_type == "pvalue") 1 else pval_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize     = min_gs_size,
    maxGSSize     = max_gs_size,
    ...
  )
  res <- filter_enrich_by_pval(res, pval_type = pval_type, pval_cutoff = pval_cutoff)
  if (is.null(res)) return(list(table = data.frame(), object = NULL))
  list(table = as.data.frame(res), object = res)
}

#' KEGG 基因集富集分析（GSEA，基于 clusterProfiler::gseKEGG）
#' @param ranked 全基因排序向量（名字为 symbol，按 log2FC 降序）
#' @param pval_type 筛选指标："padj"（默认）或 "pvalue"
gsea_kegg <- function(ranked, organism = "hsa", pval_cutoff = 0.05,
                      pval_type = c("padj", "pvalue"),
                      p_adjust_method = "BH", min_gs_size = 10, max_gs_size = 500,
                      use_internal_data = FALSE, ...) {
  pval_type     <- match.arg(pval_type)
  org_db        <- get_orgdb(organism)
  ranked_entrez <- rank_to_entrez(ranked, org_db = org_db)

  if (is.null(getOption("clusterProfiler.download.method"))) {
    options(clusterProfiler.download.method = "curl")
  }

  res <- clusterProfiler::gseKEGG(
    geneList          = ranked_entrez,
    organism          = organism,
    keyType           = "ncbi-geneid",
    pvalueCutoff      = if (pval_type == "pvalue") 1 else pval_cutoff,
    pAdjustMethod     = p_adjust_method,
    minGSSize         = min_gs_size,
    maxGSSize         = max_gs_size,
    use_internal_data = use_internal_data,
    ...
  )
  res <- filter_enrich_by_pval(res, pval_type = pval_type, pval_cutoff = pval_cutoff)
  finalize_kegg_result(res, org_db = org_db)
}

# -----------------------------------------------------------------------------
# 4. 富集结果绘图
# -----------------------------------------------------------------------------

# 内部：把 enrichResult(ORA) / gseaResult(GSEA) 统一转成绘图数据框
enrich_bubble_df <- function(enrich_obj,
                             n         = 10,
                             ont       = c("ALL", "BP", "CC", "MF"),
                             pval_type = c("padj", "pvalue")) {
  ont       <- match.arg(ont)
  pval_type <- match.arg(pval_type)
  p_col     <- if (pval_type == "padj") "p.adjust" else "pvalue"

  d <- as.data.frame(enrich_obj)
  if (nrow(d) == 0) stop("富集结果为空，无法绘图！")
  if (!p_col %in% colnames(d)) stop("富集结果缺少列：", p_col)
  d$p_val_plot <- d[[p_col]]

  if (all(c("GeneRatio", "Count") %in% colnames(d))) {
    # ORA：GeneRatio 形如 "37/301"，拆成数值；Count 直接可用
    ratio <- vapply(strsplit(as.character(d$GeneRatio), "/"),
                    function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    d$GeneRatio_num <- ratio
    d$Count_num     <- d$Count
  } else if (all(c("setSize", "core_enrichment") %in% colnames(d))) {
    # GSEA：Count = 核心富集基因数，GeneRatio = Count / setSize
    core <- as.character(d$core_enrichment)
    cnt  <- ifelse(is.na(d$core_enrichment) | core == "",
                   0, lengths(strsplit(core, "/")))
    d$Count_num     <- cnt
    d$GeneRatio_num <- cnt / d$setSize
  } else {
    stop("无法识别的富集结果：缺少 GeneRatio/Count 或 setSize/core_enrichment 列")
  }

  # 按 ont 参数筛选本体（仅当结果含 ONTOLOGY 列，即 enrich_go(ont = "ALL") 时）
  has_ont <- "ONTOLOGY" %in% colnames(d)
  if (has_ont && ont != "ALL") {
    d <- d[d$ONTOLOGY == ont, , drop = FALSE]
    if (nrow(d) == 0) stop(sprintf("结果中没有 %s 的通路！", ont))
  } else if (!has_ont && ont != "ALL") {
    warning("富集结果不含 ONTOLOGY 列（仅跑了单一本体），已忽略 ont 参数。")
  }

  # 取 top n：ont = ALL 时按本体分组各取 top n，否则整体取 top n
  if (has_ont && ont == "ALL") {
    d$ONTOLOGY <- factor(d$ONTOLOGY, levels = c("BP", "CC", "MF"))
    parts <- lapply(split(d, d$ONTOLOGY), function(g) {
      g <- g[order(g$p_val_plot, na.last = TRUE), , drop = FALSE]
      head(g, n)
    })
    d <- do.call(rbind, parts)
    # y 轴顺序：每个本体内最显著的在最上面
    lvls <- unique(unlist(lapply(parts, function(g) rev(g$Description)), use.names = FALSE))
    d$Description <- factor(d$Description, levels = lvls)
  } else {
    d <- d[order(d$p_val_plot, na.last = TRUE), , drop = FALSE]
    d <- head(d, n)
    d$Description <- factor(d$Description, levels = rev(unique(d$Description)))
  }
  d
}

#' 气泡图（ggplot2 手绘版，ORA / GSEA 通用，沿用火山图风格）
#' @param enrich_obj 富集对象
#' @param n 展示通路条数
#' @param title 图标题
#' @param ont 选择展示的 GO 本体："ALL" / "BP" / "CC" / "MF"
#' @param pval_type 颜色映射指标："padj"（默认）或 "pvalue"
#' @param colors 颜色渐变向量（默认三色过渡：蓝 -> 白 -> 红，可自定义）
#' @param facet_nrow 多本体分面行数，默认 3
plot_dot <- function(enrich_obj,
                     n          = 10,
                     title      = "",
                     ont        = c("ALL", "BP", "CC", "MF"),
                     pval_type  = c("padj", "pvalue"),
                     colors     = c("#00468BB2", "white", "#ED0000B2"),
                     facet_nrow = 3) {
  pval_type <- match.arg(pval_type)
  p_label   <- if (pval_type == "padj") "p.adjust" else "pvalue"
  d <- enrich_bubble_df(enrich_obj, n = n, ont = ont, pval_type = pval_type)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = GeneRatio_num, y = Description)) +
    ggplot2::geom_point(
      ggplot2::aes(size = Count_num, fill = p_val_plot),
      shape  = 21,
      color  = "grey30",
      stroke = 0.25,
      alpha  = 0.85
    ) +
    ggplot2::scale_fill_gradientn(colors = colors) +
    ggplot2::scale_size_continuous(range = c(3, 8)) +
    ggplot2::labs(
      x     = "GeneRatio",
      y     = NULL,
      size  = "Count",
      fill  = p_label,
      title = title
    )

  # GO 的 BP/CC/MF 分面（纵向 3 行排列，保证通路名称完整展示）
  if ("ONTOLOGY" %in% colnames(d) && length(unique(d$ONTOLOGY)) > 1) {
    p <- p + ggplot2::facet_wrap(~ ONTOLOGY, nrow = facet_nrow, scales = "free_y")
  }

  if (exists("my_theme", inherits = TRUE)) {
    p <- p + my_theme
  } else {
    p <- p + ggplot2::theme_classic(base_size = 14)
  }
  p
}

#' 柱状图（ggplot2 手绘版，ORA / GSEA 通用，支持 GO 三行分面与三色过渡）
#' @param enrich_obj 富集对象
#' @param n 展示通路条数
#' @param title 图标题
#' @param ont 选择展示的 GO 本体："ALL" / "BP" / "CC" / "MF"
#' @param pval_type 颜色映射指标："padj"（默认）或 "pvalue"
#' @param colors 颜色渐变向量（默认三色过渡：蓝 -> 白 -> 红，可自定义）
#' @param facet_nrow 多本体分面行数，默认 3
plot_bar <- function(enrich_obj,
                     n          = 10,
                     title      = "",
                     ont        = c("ALL", "BP", "CC", "MF"),
                     pval_type  = c("padj", "pvalue"),
                     colors     = c("#00468BB2", "white", "#ED0000B2"),
                     facet_nrow = 3) {
  pval_type <- match.arg(pval_type)
  p_label   <- if (pval_type == "padj") "p.adjust" else "pvalue"
  d <- enrich_bubble_df(enrich_obj, n = n, ont = ont, pval_type = pval_type)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = Count_num, y = Description, fill = p_val_plot)) +
    ggplot2::geom_col(width = 0.72, color = "grey30", linewidth = 0.2, alpha = 0.9) +
    ggplot2::scale_fill_gradientn(colors = colors) +
    ggplot2::labs(
      x     = "Count",
      y     = NULL,
      fill  = p_label,
      title = title
    )

  if ("ONTOLOGY" %in% colnames(d) && length(unique(d$ONTOLOGY)) > 1) {
    p <- p + ggplot2::facet_wrap(~ ONTOLOGY, nrow = facet_nrow, scales = "free_y")
  }

  if (exists("my_theme", inherits = TRUE)) {
    p <- p + my_theme
  } else {
    p <- p + ggplot2::theme_classic(base_size = 14)
  }
  p
}

# 内部：按 ONTOLOGY 过滤富集对象（cnet 用）
subset_enrich <- function(enrich_obj, ont) {
  res <- enrich_obj@result
  if (ont != "ALL" && "ONTOLOGY" %in% colnames(res)) {
    enrich_obj@result <- res[res$ONTOLOGY == ont, , drop = FALSE]
  }
  enrich_obj
}

#' 富集结果绘图（ORA / GSEA 通用）
#' @param enrich_obj enrich_go / enrich_kegg / gsea_go / gsea_kegg 返回的 object
#' @param p_type 绘图类型："dot" 气泡图 / "bar" 柱状图 / "cnet" 通路-基因网络图
#' @param n 展示的通路条数，默认 10
#' @param title 图标题
#' @param ont 选择展示的 GO 本体："ALL" / "BP" / "CC" / "MF"（需 enrich_go(ont = "ALL") 结果）
#' @param pval_type 显著性指标："padj"（默认）或 "pvalue"
#' @param colors 三色渐变向量，默认 c("#00468BB2", "white", "#ED0000B2")
#' @param facet_nrow 多本体分面行数，默认 3
plot_enrich <- function(enrich_obj,
                        p_type     = c("dot", "bar", "cnet"),
                        n          = 10,
                        title      = "",
                        ont        = c("ALL", "BP", "CC", "MF"),
                        pval_type  = c("padj", "pvalue"),
                        colors     = c("#00468BB2", "white", "#ED0000B2"),
                        facet_nrow = 3) {
  requireNamespace("enrichplot", quietly = TRUE)
  p_type    <- match.arg(p_type)
  ont       <- match.arg(ont)
  pval_type <- match.arg(pval_type)
  p <- switch(p_type,
              dot  = plot_dot(enrich_obj, n = n, title = title, ont = ont,
                              pval_type = pval_type, colors = colors, facet_nrow = facet_nrow),
              bar  = plot_bar(enrich_obj, n = n, title = title, ont = ont,
                              pval_type = pval_type, colors = colors, facet_nrow = facet_nrow),
              cnet = enrichplot::cnetplot(subset_enrich(enrich_obj, ont), showCategory = n))
  p
}
