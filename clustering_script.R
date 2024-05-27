# clustering of deconvoltuin data
# using heatmap and hierarchical clustering,
# the cell mean fractions are used
cell_fracs <- read.csv("cell_deconvolution_fullDF.csv")
row.names(cell_fracs) <- cell_fracs[,1]
cell_fracs <- cell_fracs[, -c(1)]

cell_fracs_mat <- as.matrix(cell_fracs, nrow = 1223, ncol = 9)


library(cluster)

dist_vals <- dist(cell_fracs_mat)
clust <- hclust(dist_vals)
plot(clust)
cut_avg <- cutree(clust, k = 9)
rect.hclust(clust , k = 9, border = 2:6)
library(dendextend)
avg_dend_obj <- as.dendrogram(clust)
avg_col_dend <- color_branches(avg_dend_obj, h = 0.8)
plot(avg_col_dend)

groups <- cutree(clust, k=9)

table(groups)
