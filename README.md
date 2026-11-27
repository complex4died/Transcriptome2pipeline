# Transcriptome Analysis Toolkit (Rebuild)

> 转录组分析工具箱（重构版）：R 语言模块化实现
>
> 核心流程：**数据清洗 → 格式转换 → 差异分析（DEG）→ 可视化（PCA / 火山图 / 热图）→ 富集分析（GO / KEGG × ORA / GSEA）→ WGCNA → 机器学习**

---

## 目录结构

```text
Rebuild/
├── R/
│   ├── functions/
│   │   ├── gene_find.R           # 基因注释查询（Ensembl ID ↔ Symbol，Human / Mouse / Rat）
│   │   ├── counts_convert.R      # 表达矩阵标准化转换（CPM / log2CPM / TPM / log2TPM）
│   │   ├── counts_degs.R         # 差异表达分析（limma-voom / DESeq2 自适应双引擎）
│   │   ├── function_plot.R       # 统一主题 my_theme + 火山图
│   │   ├── function_pca.R        # PCA 主成分分析散点图
│   │   ├── function_heatmap.R    # 差异基因热图（ComplexHeatmap 封装，支持分组切分）
│   │   ├── function_enrichment.R # 富集分析（GO / KEGG × ORA / GSEA + 绘图）
│   │   └── function_wgcna.R      # WGCNA 加权共表达网络 + 模块-性状相关性
│   ├── datas/                    # 原始输入数据（不上传 Git）
│   ├── dir_temp/                 # 测试 / 工作脚本（不上传 Git）
│   └── main.R                    # 主入口（待完善）
├── Results/                      # 分析结果输出（不上传 Git）
│   ├── Degs/
│   ├── Enrichment/
│   ├── WGCNA/
│   └── ML/
├── .gitignore
└── README.md
```

---

## 分析模块说明

### 1. `gene_find.R` — 基因信息查询

- 支持输入类型：`ENSEMBL` 或 `SYMBOL`
- 支持物种：`human`、`mouse`、`rat`（`auto` 可自动识别 Ensembl 前缀）
- 离线引擎（默认）：`org.Hs.eg.db` / `org.Mm.eg.db` / `org.Rn.eg.db`，速度毫秒级
- 在线引擎（备选）：`biomaRt`，内置容灾切换

```r
source("R/functions/gene_find.R")

info <- query_gene_info(
  genes     = exp_count$ensembl_gene_id,
  from_type = "ENSEMBL",
  species   = "auto",
  engine    = "offline"
)
```

输出列：`input_id`、`symbol`、`ensembl_id`、`gene_type`、`description`、`species`

---

### 2. `counts_convert.R` — 表达矩阵格式转换

- 输入：原始整数 Counts 矩阵（行名为基因，列名为样本）
- 支持目标格式：`"cpm"`、`"log2cpm"`（严格非负 `log2(CPM + 1)`）、`"tpm"`、`"log2tpm"`（TPM 需传 `gene_length`）

```r
source("R/functions/counts_convert.R")

# 用于 PCA / 热图 / 下游建模的标准化矩阵（严格非负，min >= 0）
exp_cpm_log1p <- convert_counts(count_matrix = count_mat, target = "log2cpm")
```

> **注意**：差异分析（DESeq2 / limma-voom）输入必须是 **Raw Counts 整数矩阵**，不能使用 CPM/TPM；而 PCA 与热图需要使用 **log 转换后的矩阵**（如 `exp_cpm_log1p`）。

---

### 3. `counts_degs.R` — 差异表达分析

- 自动适配算法：`N ≥ 30` 自动选用 `limma-voom`，`N < 30` 自动选用 `DESeq2`
- 支持多因素协变量校正（连续变量 Age + 分类变量 Gender）
- 支持 DESeq2 多核并行加速（`n_cores` 参数）
- 支持 `padj` / `pvalue` 自由切换作为筛选标准
- 阈值可自由调整：`logfc_cutoff`（默认 1.5）、`pval_cutoff`（默认 0.05）

```r
source("R/functions/counts_degs.R")

deg_out <- run_deg(
  count_matrix = count_mat,
  sample_info  = sample_info,
  group_col    = "Groups",
  case         = "sepsis",
  control      = "healthy",
  covariates   = c("Age", "Gender"),
  logfc_cutoff = 1.5,
  pval_cutoff  = 0.05,
  pval_type    = "padj",        # 或 "pvalue"
  method       = "auto",        # 或 "limma-voom" / "deseq2"
  n_cores      = 4
)
```

**返回值（list）：**

| 元素 | 类型 | 说明 |
| :--- | :--- | :--- |
| `deg_out$table` | `data.frame` | 完整差异基因大表，含 `Change` 列（Up / Down / Stable），按 `abs(log2FC)` 降序排列 |
| `deg_out$summary` | 命名向量 | `c(Total=, Up=, Down=, Stable=)`，方便取数字用于作图标注 |
| `deg_out$sig_genes` | 字符向量 | 仅显著差异基因名，直接传入下游富集分析 |

**标准化输出列名：**

| 列名 | 含义 |
| :--- | :--- |
| `symbol` | 基因名 |
| `log2FC` | Log2 Fold Change（Case vs Control） |
| `AveExpr` | 平均表达量 |
| `pvalue` | 原始 P 值 |
| `padj` | BH 校正后 P 值（FDR） |
| `Change` | 分类标注：`Up` / `Down` / `Stable` |

---

### 4. `function_plot.R` — 统一主题 + 火山图

- `my_theme`：全局 ggplot2 主题（白边半透明圆点、零点虚线、黑边面板）
- `plot_volcano()`：支持 `padj` / `pvalue` 切换、自定义阈值虚线、对称 X 轴、Top 基因标注（`ggrepel`）

```r
source("R/functions/function_plot.R")

p_volcano <- plot_volcano(
  deg_data     = deg_out$table,
  pval_type    = "padj",
  pval_cutoff  = 0.05,
  logfc_cutoff = 1.5,
  top_n        = 10,
  title        = "Volcano Plot (sepsis vs healthy)"
)
p_volcano
```

---

### 5. `function_pca.R` — PCA 主成分分析散点图

- 沿用火山图统一视觉风格（`shape = 21` 白边半透明圆点、零点虚线参考线、`my_theme` 主题）
- 自动对齐 `exp_mat` 与 `sample_info$Sample_ID`，并过滤零方差 / `NA` 基因
- 坐标轴自动标注 `PC1`、`PC2` 方差解释率
- 默认配色直接暴露在函数签名中：`group_colors = c("control" = "#00468BB2", "case" = "#ED0000B2")`
- 支持 `show_ellipse = TRUE`（95% 置信椭圆）与 `top_n_var`（仅用高变基因计算 PCA）

```r
source("R/functions/function_pca.R")

p_pca <- plot_pca(
  exp_mat      = exp_cpm_log1p,
  sample_info  = sample_info,
  group_col    = "Groups",
  case         = "sepsis",
  control      = "healthy",
  top_n_var    = NULL,          # 可选，如 5000 表示取方差前 5000 个基因
  show_ellipse = TRUE,
  title        = "PCA: Sepsis vs Healthy"
)
p_pca
```

---

### 6. `function_heatmap.R` — 差异基因表达热图

- 基于 `ComplexHeatmap` 封装，自动提取上下调各 `top_n` 个显著差异基因（默认 `top_n = 50`，共 100 个基因）
- 自动按行做 Z-score 标准化（`-2` 蓝 → `0` 白 → `2` 红）
- `split_by_group = TRUE`（默认）：按 `control`（左）与 `case`（右）分组切分列并保留组内聚类；设为 `FALSE` 可恢复全样本无监督聚类
- 自动对齐 `sample_info$Sample_ID` 并动态绑定 `control`/`case` 分组颜色

```r
source("R/functions/function_heatmap.R")

p_heat <- plot_heatmap(
  exp_mat        = exp_cpm_log1p,
  deg_table      = deg_out$table,
  sample_info    = sample_info,
  case           = "sepsis",
  control        = "healthy",
  top_n          = 25,
  split_by_group = TRUE,
  title          = "Top DEGs: Sepsis vs Healthy"
)
ComplexHeatmap::draw(p_heat)
```

---

### 7. `function_enrichment.R` — 富集分析（GO / KEGG × ORA / GSEA）

统一封装 **过表达分析（ORA）** 与 **基因集富集分析（GSEA）**，GO / KEGG 各一套；输入输出均为 gene symbol。

**输入准备函数：**

| 函数 | 作用 |
| :--- | :--- |
| `get_sig_genes(deg_table, direction)` | 提取显著基因（`both` / `up` / `down`），ORA 输入 |
| `get_universe(exp_mat)` | 提取背景基因集（全部基因），ORA universe |
| `rank_genes(deg_table)` | 生成按 log2FC 降序的全基因排序向量，GSEA 输入 |

**分析函数（`pval_type` 支持 `padj` / `pvalue` 全局切换）：**

| 函数 | 类型 | 说明 |
| :--- | :--- | :--- |
| `enrich_go(genes, universe, ont)` | ORA | GO 过表达，`ont` 支持 `BP`/`MF`/`CC`/`ALL` |
| `enrich_kegg(genes, universe)` | ORA | KEGG 过表达 |
| `gsea_go(ranked, ont)` | GSEA | GO 基因集富集 |
| `gsea_kegg(ranked)` | GSEA | KEGG 基因集富集 |

**绘图函数：**

| 函数 | 说明 |
| :--- | :--- |
| `plot_dot(enrich_obj, n, ...)` | 气泡图（ORA / GSEA 通用） |
| `plot_bar(enrich_obj, n, ...)` | 柱状图 |
| `plot_enrich(enrich_obj, p_type, ...)` | 统一入口：`dot` / `bar` / `cnet` |

```r
source("R/functions/function_enrichment.R")

sig_genes <- get_sig_genes(deg_out$table, direction = "both")
universe  <- get_universe(count_mat)
ranked    <- rank_genes(deg_out$table)

go_ora  <- enrich_go(sig_genes, universe, ont = "ALL", pval_type = "padj")
gsea_go <- gsea_go(ranked, pval_type = "padj")

if (nrow(go_ora$table) > 0) {
  print(plot_enrich(go_ora$object, p_type = "dot", n = 5, ont = "ALL", title = "GO ORA"))
}
```

> **注意**：`enrich_kegg` / `gsea_kegg` 依赖 KEGG 在线 REST API；若网络不通，会返回空结果或报错，此时可用 GO 结果或改用本地 MSigDB 基因集。

---

### 8. `function_wgcna.R` — WGCNA 加权共表达网络

封装 **网络构建 → 模块识别 → 模块-性状相关性** 全流程，关键参数全部暴露在签名中。

```r
source("R/functions/function_wgcna.R")

wgcna_res <- run_wgcna(
  exp_mat          = exp_cpm_log1p,   # log 处理后的矩阵
  sample_info      = sample_info,
  group_col        = "Groups",
  case             = "sepsis",
  control          = "healthy",
  trait_cols       = c("Age", "Gender"),
  n_var_genes      = 8000,            # 保留方差最大的前 N 基因
  network_type     = "signed",
  soft_power       = NULL,            # NULL = 按无标度拓扑自动选择
  n_threads        = 4,               # 并行线程数
  min_module_size  = 30,
  merge_cut_height = 0.25
)

print(table(wgcna_res$module_colors))
plot_module_trait(wgcna_res, title = "WGCNA Module-Trait Relationships")
```

**返回值（list）**：`dat_expr`、`soft_power`、`net`、`module_colors`、`MEs`、`traits`、`module_trait_cor`、`module_trait_pval`

---

## ⚠️ 编码规范与已知注意事项

### 1. `select()` 函数命名冲突（Bioconductor vs dplyr）

加载 `org.Hs.eg.db`、`AnnotationDbi` 或 `clusterProfiler` 后，  
Bioconductor 的 `AnnotationDbi::select` 会覆盖 `dplyr::select`，  
导致管道中的 `select()` 报错：
```text
unable to find an inherited method for function 'select' for signature '"data.frame"'
```

**全项目统一规范**：涉及列筛选时，始终显式写为 `dplyr::select(...)`。

```r
# ❌ 易冲突
df %>% select(symbol, log2FC)

# ✅ 显式指定
df %>% dplyr::select(symbol, log2FC)
```

### 2. `rename()` 函数命名冲突（S4Vectors vs dplyr）

加载 `DESeq2`（依赖 `S4Vectors`）后，`S4Vectors::rename` 会覆盖 `dplyr::rename`，  
导致 `rename(log2FC = logFC)` 报错 `object 'logFC' not found`。  
**规范**：列重命名统一写为 `dplyr::rename(...)`（已在 `counts_degs.R` 中修正）。

### 3. `cor()` S4 泛型冲突（S4Vectors/IRanges vs WGCNA）

`S4Vectors`/`IRanges` 会把 `stats::cor` 变成 S4 泛型；WGCNA 内部通过 `match.fun("cor")` 查找 `cor`，  
若不 `attach` WGCNA 就会命中该泛型并报错 `unused arguments (weights.x, weights.y, cosine)`。  
**规范**：`run_wgcna()` 内部已用 `library(WGCNA)` 显式 attach，确保 `cor` 命中 WGCNA 版本。

### 4. 差异分析算法选择依据（样本量标准）

| 场景 | 推荐算法 | 原因 |
| :--- | :--- | :--- |
| `N ≥ 30`（大样本临床队列） | **`limma-voom`** | 闭式解，400 样本只需约 6 秒；天然支持连续协变量 |
| `N < 30`（小样本实验设计） | **`DESeq2`** | 经验贝叶斯方差收缩，小样本下统计更可靠 |

### 5. DESeq2 / edgeR 输入必须是原始整数 Counts

`DESeq2` 和 `edgeR` 底层基于负二项分布回归，  
严禁将 CPM / TPM / log-transformed 矩阵直接传入，否则结果错误。

### 6. `ggplot2` 统一主题（`my_theme`）配色与坐标轴规范

- 全局配色基调：Control / Down `#00468BB2`，Stable `grey60`，Case / Up `#ED0000B2`
- `my_theme` 使用 `panel.border` 绘制黑边面板，网格关闭，标题居中加粗

---

## 依赖包清单

### CRAN 包
```r
install.packages(c("dplyr", "ggplot2", "ggsci", "ggrepel",
                   "circlize", "data.table", "tibble"))
```

### Bioconductor 包
```r
BiocManager::install(c(
  "limma", "edgeR", "DESeq2", "ComplexHeatmap",
  "clusterProfiler", "enrichplot", "AnnotationDbi",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
  "WGCNA", "BiocParallel", "biomaRt"
))
```

---

## 数据输入规范

| 文件 | 格式要求 |
| :--- | :--- |
| `exp_count.csv` | 行名为 Gene Symbol，列名为样本 ID，值为非负整数 |
| `sample_info.csv` | 必须包含 `Sample_ID`、`Groups` 列；可含 `Age`、`Gender` 等协变量列 |

---

## 快速上手

完整流程见 `R/dir_temp/test.R`，核心步骤：

```r
rm(list = ls())

# 一键加载所有功能模块
func_dir <- "R/functions/"
invisible(lapply(list.files(func_dir, pattern = "\\.R$", full.names = TRUE), source))

# 读取数据
library(data.table)
library(tibble)
exp_count   <- fread("R/datas/exp_count.csv") %>% column_to_rownames("symbol")
sample_info <- fread("R/datas/sample_info.csv")
count_mat   <- as.matrix(exp_count)

# 1. 格式转换（用于 PCA / 热图 / WGCNA）
exp_cpm_log1p <- convert_counts(count_matrix = count_mat, target = "log2cpm")

# 2. 差异分析
deg_out <- run_deg(count_mat, sample_info,
                   case = "sepsis", control = "healthy",
                   covariates = c("Age", "Gender"))
deg_out$summary

# 3. 可视化：火山图 + PCA + 热图
p_volcano <- plot_volcano(deg_out$table, title = "Volcano Plot")
p_pca     <- plot_pca(exp_cpm_log1p, sample_info, case = "sepsis", control = "healthy")
p_heat    <- plot_heatmap(exp_cpm_log1p, deg_out$table, sample_info,
                          case = "sepsis", control = "healthy", top_n = 25)

# 4. 富集分析
sig_genes <- get_sig_genes(deg_out$table, direction = "both")
universe  <- get_universe(count_mat)
ranked    <- rank_genes(deg_out$table)
go_ora    <- enrich_go(sig_genes, universe, ont = "ALL", pval_type = "padj")
gsea_go   <- gsea_go(ranked, pval_type = "padj")

# 5. WGCNA
wgcna_res <- run_wgcna(exp_cpm_log1p, sample_info,
                       case = "sepsis", control = "healthy",
                       trait_cols = c("Age", "Gender"))
plot_module_trait(wgcna_res)
```
