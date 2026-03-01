set.seed(42)
genes <- c(paste0("TF", 1:8), paste0("Gene", 1:42))
cells <- paste0("Cell", 1:120)
expr <- matrix(rnorm(length(genes) * length(cells), mean = 0, sd = 1),
               nrow = length(genes),
               dimnames = list(genes, cells))
# add subgroup signal
expr[1:5, 1:60] <- expr[1:5, 1:60] + 0.8
expr[6:10, 61:120] <- expr[6:10, 61:120] + 0.8

singlecell_data <- expr
save(singlecell_data, file = "data/singlecell_data.rda", compress = "xz")
