library(dplyr)

#' 多物种通用基因信息查询工具（支持 人 / 小鼠 / 大鼠）
#'
#' @param genes 字符向量，输入的基因列表
#' @param from_type 输入的 ID 类型："ENSEMBL" 或 "SYMBOL"
#' @param species 物种："auto"（自动推断）、"human"（人）、"mouse"（小鼠）、"rat"（大鼠）
#' @param engine 注释引擎："offline"（推荐）或 "online"
#' @return 标准化数据框 (input_id, symbol, ensembl_id, gene_type, description, species)
query_gene_info <- function(genes, 
                            from_type = c("ENSEMBL", "SYMBOL"), 
                            species = c("auto", "human", "mouse", "rat"), 
                            engine = "offline") {
  from_type <- match.arg(from_type)
  species   <- match.arg(species)
  
  # 1. 基础清洗：去除空格与 Ensembl 版本号
  genes_clean <- trimws(as.character(genes))
  if (from_type == "ENSEMBL") {
    genes_clean <- gsub("\\..*", "", genes_clean)
  }
  genes_unique <- unique(genes_clean[genes_clean != "" & !is.na(genes_clean)])
  if (length(genes_unique) == 0) stop("输入的有效基因列表为空！")
  # 2. 智能物种识别
  if (species == "auto") {
    if (from_type == "ENSEMBL") {
      sample_id <- genes_unique[1]
      if (grepl("^ENSG", sample_id)) {
        species <- "human"
      } else if (grepl("^ENSMUSG", sample_id)) {
        species <- "mouse"
      } else if (grepl("^ENSRNOG", sample_id)) {
        species <- "rat"
      } else {
        species <- "human"
      }
      message(sprintf(">>> [Auto-Detect] 检测到 Ensembl 前缀，已自动判定物种为: 【%s】", species))
    } else {
      species <- "human"
    }
  }
  # 3. 物种配置映射
  spec_config <- list(
    human = list(pkg = "org.Hs.eg.db", mart_dataset = "hsapiens_gene_ensembl", symbol_attr = "hgnc_symbol"),
    mouse = list(pkg = "org.Mm.eg.db", mart_dataset = "mmusculus_gene_ensembl", symbol_attr = "mgi_symbol"),
    rat   = list(pkg = "org.Rn.eg.db", mart_dataset = "rnorvegicus_gene_ensembl", symbol_attr = "rgd_symbol")
  )
  cfg <- spec_config[[species]]
  message(sprintf(">>> [GeneQuery] 正在查询 %d 个基因 (物种: %s, 类型: %s, 引擎: %s)...", 
                  length(genes_unique), species, from_type, engine))
  # =========================================================================
  # 引擎 1：本地离线引擎（已修复重命名冲突 Bug）
  # =========================================================================
  if (engine == "offline") {
    if (!requireNamespace(cfg$pkg, quietly = TRUE)) {
      stop(sprintf("请先安装对应物种的注释包：BiocManager::install('%s')", cfg$pkg))
    }
    suppressPackageStartupMessages(library(cfg$pkg, character.only = TRUE))
    db_obj <- get(cfg$pkg)
    
    raw_res <- AnnotationDbi::select(
      db_obj,
      keys     = genes_unique,
      columns  = c("ENSEMBL", "SYMBOL", "GENETYPE", "GENENAME"),
      keytype  = from_type
    )
    
    # 修复关键：用 mutate 显式赋值，避免 rename 同名列冲突
    res <- raw_res %>%
      mutate(
        input_id    = .data[[from_type]],
        ensembl_id  = ENSEMBL,
        symbol      = SYMBOL,
        gene_type   = GENETYPE,
        description = GENENAME,
        species     = species
      ) %>%
      dplyr::select(input_id, symbol, ensembl_id, gene_type, description, species) %>%
      distinct(input_id, .keep_all = TRUE)
    
    message(">>> 查询完成！")
    return(res)
  }
  # =========================================================================
  # 引擎 2：在线引擎 (biomaRt)
  # =========================================================================
  if (engine == "online") {
    if (!requireNamespace("biomaRt", quietly = TRUE)) stop("请先安装 biomaRt")
    library(biomaRt)
    
    hosts <- c("https://www.ensembl.org", "https://dec2021.archive.ensembl.org")
    mart <- NULL
    for (h in hosts) {
      mart <- tryCatch({
        useMart("ensembl", dataset = cfg$mart_dataset, host = h)
      }, error = function(e) NULL)
      if (!is.null(mart)) break
    }
    if (is.null(mart)) stop("Ensembl 在线节点响应超时，请切换为 engine = 'offline'")
    
    filter_key <- if (from_type == "ENSEMBL") "ensembl_gene_id" else cfg$symbol_attr
    
    raw_res <- getBM(
      attributes = c("ensembl_gene_id", cfg$symbol_attr, "gene_biotype", "description"),
      filters    = filter_key,
      values     = genes_unique,
      mart       = mart
    )
    
    colnames(raw_res)[colnames(raw_res) == cfg$symbol_attr] <- "symbol"
    res <- raw_res %>%
      rename(
        ensembl_id  = ensembl_gene_id,
        gene_type   = gene_biotype,
        description = description
      ) %>%
      mutate(
        input_id = if (from_type == "ENSEMBL") ensembl_id else symbol,
        species  = species
      ) %>%
      dplyr::select(input_id, symbol, ensembl_id, gene_type, description, species) %>%
      distinct(input_id, .keep_all = TRUE)
    
    message(">>> 查询完成！")
    return(res)
  }
}