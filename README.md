# GRNSimilarity

用于单细胞基因调控网络（GRN）分析的 R 包，支持：

- 使用 **hdWGCNA** / **corto(ARACNE)** / **SCENIC** 推断网络；
- 进行转录因子（TF）虚拟敲除（virtual KO）；
- 计算敲除前后网络相似度变化（支持 `cosine` / `pearson` / `jaccard` / `spearman` / `frobenius`）；
- 在两个细胞亚群间对比 TF 对网络扰动的影响并排序；
- 输出可直接用于 `ggplot2` 的长表数据（`plot_data`）。

## 依赖包与方法参数

包已集成以下后端依赖（可选）：

- `hdWGCNA`（方法 `method = "hdWGCNA"`）
- `corto`（方法 `method = "ARACNE"`）
- `SCENIC`（方法 `method = "SCENIC"`）

通过 `inferGRN(..., use_external = TRUE, method_params = list(...))` 设置参数：

- `hdwgcna_network_type`: `"unsigned"` / `"signed"`
- `corto_dpi_tolerance`: 默认 `0.1`
- `corto_nboot`: 默认 `100`
- `scenic_corr_threshold`: 默认 `0.03`

## 核心函数

- `inferGRN(expr_data, method, power, k, regulon_db, use_external, method_params)`
- `simulateKO(grn, tf_list)`
- `calculateSimilarity(grn_before, grn_after, metrics, edge_threshold)`
- `networkRanking(expr_data, cells_col, cells_A, cells_B, tf_list, method, similarity_metrics, primary_metric)`

## 使用示例

```r
library(GRNSimilarity)

expr <- singlecell_data()
cells_col <- c(rep("A", 60), rep("B", 60))
tf_list <- c("TF1", "TF2", "TF3")

ranking_results <- networkRanking(
  expr_data = expr,
  cells_col = cells_col,
  cells_A = "A",
  cells_B = "B",
  tf_list = tf_list,
  method = "hdWGCNA",
  use_external = TRUE,
  method_params = list(
    hdwgcna_network_type = "unsigned",
    corto_dpi_tolerance = 0.1,
    corto_nboot = 100,
    scenic_corr_threshold = 0.03
  ),
  similarity_metrics = c("cosine", "jaccard", "pearson"),
  primary_metric = "cosine"
)

# 1) TF 排名（主指标）
print(ranking_results$ranking)

# 2) 所有指标的长表（适合 ggplot2）
head(ranking_results$plot_data$tf_ranking)

# 3) 基因分数（适合条形图）
head(ranking_results$plot_data$gene_scores)
```
