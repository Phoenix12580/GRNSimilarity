# GRNSimilarity

用于单细胞基因调控网络（GRN）分析的 R 包，支持：

- 使用 **hdWGCNA** / **ARACNE** / **SCENIC** 风格方法推断网络；
- 进行转录因子（TF）虚拟敲除（virtual KO）；
- 计算敲除前后网络相似度变化；
- 在两个细胞亚群间对比 TF 对网络扰动的影响并排序。

## 核心函数

- `inferGRN(expr_data, method, power, k, regulon_db)`
- `simulateKO(grn, tf_list)`
- `calculateSimilarity(grn_before, grn_after)`
- `networkRanking(expr_data, cells_col, cells_A, cells_B, tf_list, method)`

## 使用示例

```r
library(GRNSimilarity)

expr <- singlecell_data()
cells_col <- c(rep("A", 60), rep("B", 60))
tf_list <- c("TF1", "TF2")

ranking_results <- networkRanking(
  expr_data = expr,
  cells_col = cells_col,
  cells_A = "A",
  cells_B = "B",
  tf_list = tf_list,
  method = "hdWGCNA"
)

print(ranking_results$ranking)
```

