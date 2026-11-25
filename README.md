# Transcriptome Analysis Toolkit (Rebuild)

> 转录组分析工具箱（重构版）：R 语言模块化实现
> 
> 核心流程：**数据清洗 → 格式转换 → 差异分析（DEG）→ 富集分析 → WGCNA → 机器学习**

---

## 目录结构

```text
Rebuild/
├── R/
│   └── functions/
│       ├── gene_find.R          # 基因注释查询（Ensembl ID ↔ Symbol，支持 Human / Mouse / Rat）
│       ├── counts_convert.R     # 表达矩阵标准化转换（CPM / log2CPM / TPM / log2TPM）
│       ├── deg.R                # 差异表达分析（limma-voom / DESeq2 自适应双引擎）
│       ├── enrichment.R         # 功能富集分析（GO / KEGG / GSEA）
│       ├── wgcna.R              # 加权基因共表达网络分析
│       └── ml.R                 # 机器学习特征筛选（LASSO + RF + SVM-RFE）
├── datas/                       # 原始输入数据（不上传 Git）
├── Results/                     # 分析结果输出（不上传 Git）
│   ├── Degs/
│   ├── Enrichment/
│   ├── WGCNA/
│   └── ML/
├── .gitignore
└── README.md
```

---

## 分析模块说明

### `gene_find.R` — 基因信息查询

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

### `counts_convert.R` — 表达矩阵格式转换

- 输入：原始整数 Counts 矩阵（行名为基因，列名为样本）
- 支持目标格式：`"cpm"`、`"log2cpm"`（严格非负）、`"tpm"`、`"log2tpm"`

```r
source("R/functions/counts_convert.R")

# 用于绘图的标准化矩阵（非负，min=0）
exp_plot <- convert_counts(count_matrix = exp_mat, target = "log2cpm")
```

> **注意**：差异分析（DESeq2 / limma-voom）输入必须是 **Raw Counts 整数矩阵**，不能使用 CPM/TPM！

---

### `deg.R` — 差异表达分析

- 自动适配算法：`N ≥ 30` 自动选用 `limma-voom`，`N < 30` 自动选用 `DESeq2`
- 支持多因素协变量校正（连续变量 Age + 分类变量 Gender）
- 支持 DESeq2 多核并行加速（`n_cores` 参数）
- 支持 `padj` / `pvalue` 自由切换作为筛选标准
- 阈值可自由调整：`logfc_cutoff`（默认 1.5）、`pval_cutoff`（默认 0.05）

```r
source("R/functions/deg.R")

deg_out <- run_deg(
  count_matrix = exp_mat,
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

### 2. 差异分析算法选择依据（样本量标准）

| 场景 | 推荐算法 | 原因 |
| :--- | :--- | :--- |
| `N ≥ 30`（大样本临床队列） | **`limma-voom`** | 闭式解，400样本只需约 6 秒；天然支持连续协变量 |
| `N < 30`（小样本实验设计） | **`DESeq2`** | 经验贝叶斯方差收缩，小样本下统计更可靠 |

### 3. DESeq2 / edgeR 输入必须是原始整数 Counts

`DESeq2` 和 `edgeR` 底层基于负二项分布回归，  
严禁将 CPM / TPM / log-transformed 矩阵直接传入，否则结果错误。

---

## 依赖包清单

### CRAN 包
```r
install.packages(c("dplyr", "ggplot2", "ggsci", "ggrepel", 
                   "ComplexHeatmap", "circlize", "data.table", "tibble"))
```

### Bioconductor 包
```r
BiocManager::install(c(
  "limma", "edgeR", "DESeq2",
  "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
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

```r
rm(list = ls())

# 一键加载所有功能模块
path_function <- "R/functions/"
invisible(lapply(list.files(path_function, pattern = "\\.R$", full.names = TRUE), source))

# 读取数据
library(data.table)
library(tibble)
exp_count  <- fread("datas/exp_count.csv")  %>% column_to_rownames("symbol")
sample_info <- fread("datas/sample_info.csv")

# 格式转换（用于绘图）
exp_plot <- convert_counts(count_matrix = as.matrix(exp_count), target = "log2cpm")

# 差异分析
deg_out <- run_deg(
  count_matrix = as.matrix(exp_count),
  sample_info  = sample_info,
  case         = "sepsis",
  control      = "healthy",
  covariates   = c("Age", "Gender")
)

# 查看统计摘要
deg_out$summary
```
