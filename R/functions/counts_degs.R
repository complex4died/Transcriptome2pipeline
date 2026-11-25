library(edgeR)
library(limma)
library(dplyr)

#' 通用差异表达分析模块（支持 Pval/Padj 自由切换版）
#'
#' @param count_matrix 纯数值表达矩阵（行名为基因，列名为样本）
#' @param sample_info 临床信息数据框
#' @param group_col 分组列名，如 "Groups"
#' @param case 实验组/疾病组名称，如 "sepsis"
#' @param control 对照组/健康组名称，如 "healthy"
#' @param covariates 协变量向量（如 c("Age", "Gender")，NULL 为不校正）
#' @param logfc_cutoff log2FC 绝对值阈值，默认 1.5
#' @param pval_cutoff P值筛选阈值，默认 0.05
#' @param pval_type 筛选使用的统计指标："padj"（默认，FDR校正P值）或 "pvalue"（原始P值）
#' @param method "auto"（>=30走voom，<30走deseq2）、"limma-voom" 或 "deseq2"
#' @param n_cores 并行计算核心数，默认 4
#' @return list 包含 3 个元素：table, summary, sig_genes
run_deg <- function(count_matrix, 
                    sample_info, 
                    group_col    = "Groups", 
                    case         = "sepsis", 
                    control      = "healthy", 
                    covariates   = NULL, 
                    logfc_cutoff = 1.5,
                    pval_cutoff  = 0.05,
                    pval_type    = c("padj", "pvalue"),
                    method       = "auto",
                    n_cores      = 4) {
  
  pval_type <- match.arg(pval_type) # 严格校验参数为 "padj" 或 "pvalue"
  
  # -------------------------------------------------------------------------
  # 1. 样本对齐与合法性检查
  # -------------------------------------------------------------------------
  common_samples <- intersect(colnames(count_matrix), sample_info$Sample_ID)
  if (length(common_samples) == 0) {
    stop("表达矩阵列名与 sample_info 的 Sample_ID 无任何重合！请核对样本名。")
  }
  
  mat  <- count_matrix[, common_samples]
  meta <- sample_info[match(common_samples, sample_info$Sample_ID), ]
  
  # 确保基准水平（control 在前，case 在后）
  meta[[group_col]] <- factor(meta[[group_col]], levels = c(control, case))
  
  # -------------------------------------------------------------------------
  # 2. 协变量缺失值处理
  # -------------------------------------------------------------------------
  if (!is.null(covariates)) {
    valid_covs <- intersect(covariates, colnames(meta))
    if (length(valid_covs) != length(covariates)) {
      warning("部分协变量未找到，仅使用有效协变量: ", paste(valid_covs, collapse = ", "))
    }
    complete_idx <- complete.cases(meta[, c(group_col, valid_covs)])
    if (sum(!complete_idx) > 0) {
      message(sprintf(">>> 检测到协变量存在 NA，已自动剔除 %d 个样本...", sum(!complete_idx)))
      mat  <- mat[, complete_idx]
      meta <- meta[complete_idx, ]
    }
    covariates <- valid_covs
  }
  
  # -------------------------------------------------------------------------
  # 3. 算法自适应选择
  # -------------------------------------------------------------------------
  n_samples <- ncol(mat)
  if (method == "auto") {
    method <- if (n_samples >= 30) "limma-voom" else "deseq2"
  }
  message(sprintf(">>> [DEG] 样本量: %d | 算法: %s | 比较组: %s (Case) vs %s (Control)", 
                  n_samples, method, case, control))
  if (!is.null(covariates) && length(covariates) > 0) {
    message(">>> [DEG] 纳入多因素校正: ", paste(covariates, collapse = " + "))
  }
  
  raw_deg <- NULL
  
  # =========================================================================
  # 引擎 A：limma-voom (大样本模式)
  # =========================================================================
  if (method == "limma-voom") {
    dge <- DGEList(counts = mat)
    dge <- calcNormFactors(dge)
    
    if (!is.null(covariates) && length(covariates) > 0) {
      formula_str <- paste("~ 0 +", group_col, "+", paste(covariates, collapse = " + "))
    } else {
      formula_str <- paste("~ 0 +", group_col)
    }
    design <- model.matrix(as.formula(formula_str), data = meta)
    
    colnames(design) <- make.names(colnames(design))
    col_case    <- make.names(paste0(group_col, case))
    col_control <- make.names(paste0(group_col, control))
    
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    
    contrast_str <- paste0(col_case, " - ", col_control)
    contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)
    
    fit2 <- contrasts.fit(fit, contrast_mat)
    fit2 <- eBayes(fit2)
    
    raw_deg <- topTable(fit2, number = Inf, adjust.method = "BH") %>%
      tibble::rownames_to_column(var = "symbol") %>%
      rename(
        log2FC  = logFC,      # 统一 Log Fold 2 Change
        AveExpr = AveExpr,
        pvalue  = P.Value,    # 统一 P.value
        padj    = adj.P.Val   # 统一 P.adjust
      )
  }
  
  # =========================================================================
  # 引擎 B：DESeq2 (小样本模式 + 多核心并行加速)
  # =========================================================================
  if (method == "deseq2") {
    if (!requireNamespace("DESeq2", quietly = TRUE)) stop("请先安装 DESeq2")
    library(DESeq2)
    
    if (!is.null(covariates) && length(covariates) > 0) {
      design_formula <- as.formula(paste("~", paste(covariates, collapse = " + "), "+", group_col))
    } else {
      design_formula <- as.formula(paste("~", group_col))
    }
    
    dds <- DESeqDataSetFromMatrix(countData = round(mat), colData = meta, design = design_formula)
    
    use_parallel <- FALSE
    bpparam <- NULL
    if (n_cores > 1 && requireNamespace("BiocParallel", quietly = TRUE)) {
      use_parallel <- TRUE
      message(sprintf(">>> [DESeq2] 已开启多核加速，正在使用 %d 个 CPU 核心并行计算...", n_cores))
      bpparam <- if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(workers = n_cores)
      } else {
        BiocParallel::MulticoreParam(workers = n_cores)
      }
    }
    
    dds <- DESeq(dds, parallel = use_parallel, BPPARAM = bpparam)
    res <- results(dds, contrast = c(group_col, case, control), parallel = use_parallel, BPPARAM = bpparam)
    
    raw_deg <- as.data.frame(res) %>%
      tibble::rownames_to_column(var = "symbol") %>%
      rename(
        log2FC  = log2FoldChange,
        AveExpr = baseMean,
        pvalue  = pvalue,
        padj    = padj
      )
  }
  
  # =========================================================================
  # 4. 核心清洗：根据 pval_type 进行 Change 标注 + 绝对值降序排列
  # =========================================================================
  # 动态指定使用 padj 还是 pvalue 判定阈值
  stat_col <- if (pval_type == "padj") "padj" else "pvalue"
  
  deg_clean <- raw_deg %>%
    mutate(
      Change = case_when(
        !is.na(.data[[stat_col]]) & .data[[stat_col]] < pval_cutoff & log2FC >= logfc_cutoff  ~ "Up",
        !is.na(.data[[stat_col]]) & .data[[stat_col]] < pval_cutoff & log2FC <= -logfc_cutoff ~ "Down",
        TRUE ~ "Stable"
      )
    ) %>%
    # 核心：按 log2FC 绝对值从大到小严格降序排列
    arrange(desc(abs(log2FC))) %>%
    dplyr::select(symbol, log2FC, AveExpr, pvalue, padj, Change, everything())
  
  # -------------------------------------------------------------------------
  # 5. 统计与多元素打包返回
  # -------------------------------------------------------------------------
  n_up     <- sum(deg_clean$Change == "Up")
  n_down   <- sum(deg_clean$Change == "Down")
  n_stable <- sum(deg_clean$Change == "Stable")
  
  message("=================================================================")
  message(sprintf(">>> [统计摘要] 总基因数: %d | Up: %d | Down: %d | Stable: %d", 
                  nrow(deg_clean), n_up, n_down, n_stable))
  message(sprintf(">>> [当前阈值] |log2FC| >= %.2f 且 %s < %.3f", 
                  logfc_cutoff, pval_type, pval_cutoff))
  
  # 贴心提醒：如果严格的 padj 筛不出基因，自动给出提示
  if (n_up + n_down == 0 && pval_type == "padj") {
    message(">>> [提示] 当前在 padj 严格阈值下未筛选出显著基因，建议尝试传入 pval_type = 'pvalue' 进行探索！")
  }
  message("=================================================================")
  
  summary_counts <- c(
    Total  = nrow(deg_clean),
    Up     = n_up,
    Down   = n_down,
    Stable = n_stable
  )
  
  sig_genes_vec <- deg_clean %>% 
    dplyr::filter(Change %in% c("Up", "Down")) %>% 
    dplyr::pull(symbol)
  
  return(list(
    table     = deg_clean,       # 主数据框（无论选哪种指标，pvalue 和 padj 双列都完整保留）
    summary   = summary_counts,  # 命名统计向量
    sig_genes = sig_genes_vec    # 纯显著差异基因字符向量
  ))
}