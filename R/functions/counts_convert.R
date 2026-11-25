#' Count 表达矩阵标准化转换小工具
#'
#' @param count_matrix 纯数值矩阵（行名为基因，列名为样本）
#' @param target 目标格式："cpm", "log2cpm", "tpm", "log2tpm"
#' @param gene_length 命名数值向量，每个基因的长度(bp)，计算 TPM 时必须提供
#' @return 标准化后的矩阵
convert_counts <- function(count_matrix, 
                           target = c("log2cpm", "cpm", "tpm", "log2tpm"), 
                           gene_length = NULL) {
  target <- match.arg(target)
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("请先安装 edgeR：BiocManager::install('edgeR')")
  }
  
  # 1. CPM 与 log2-CPM 计算
  if (target == "cpm") {
    message(">>> [Convert] 正在计算 CPM...")
    return(edgeR::cpm(count_matrix, log = FALSE))
  }
  
  # if (target == "log2cpm") {
  #  message(">>> [Convert] 正在计算 log2(CPM + 1)...")
  #  return(edgeR::cpm(count_matrix, log = TRUE, prior.count = 1))
  # } 
  
  if (target == "log2cpm") {
    message(">>> [Convert] 正在计算严格非负的 log2(CPM + 1)...")
    cpm_mat <- edgeR::cpm(count_matrix, log = FALSE)
    return(log2(cpm_mat + 1))
  }
  
  
  # 2. TPM 与 log2-TPM 计算（需要基因长度）
  if (target %in% c("tpm", "log2tpm")) {
    if (is.null(gene_length)) {
      stop("计算 TPM 必须传入 gene_length (基因长度) 参数！若无长度数据，建议使用 target = 'log2cpm'")
    }
    
    # 确保基因长度与矩阵行名对齐
    common_genes <- intersect(rownames(count_matrix), names(gene_length))
    if (length(common_genes) == 0) stop("gene_length 的基因名与表达矩阵行名无任何重合！")
    
    mat <- count_matrix[common_genes, , drop = FALSE]
    len <- gene_length[common_genes]
    
    message(sprintf(">>> [Convert] 匹配到 %d 个基因的长度信息，正在计算 TPM...", length(common_genes)))
    
    # TPM 计算核心公式
    rate <- mat / (len / 1000)
    tpm_mat <- sweep(rate, 2, colSums(rate), "/") * 1e6
    
    if (target == "log2tpm") {
      return(log2(tpm_mat + 1))
    } else {
      return(tpm_mat)
    }
  }
}