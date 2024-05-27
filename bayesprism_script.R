## Deconvolution using BayesPrism
library(BayesPrism)

# reading data
# single cell data
scrnaseq <- readRDS("./wuetal_scseq.rds")
scseq_mat <- scrnaseq[["RNA"]]$counts
scrnaseq_mat <- t(as.matrix(scseq_mat)) # reformating to required way
saveRDS(scrnaseq_mat, "sc_matrix.rds")

# bulk data
bulk_rna <- read.csv("raw_counts.tsv", header = TRUE,
                     sep = "\t") # raw counts only, get gene symbols
bulk_rna_mat <- t(bulk_rna)
colnames(bulk_rna_mat) <- bulk_rna_mat[1,]
bulk_rna_mat <- bulk_rna_mat[2:1224,]
bulk_rna_mat2 <- matrix(as.numeric(bulk_rna_mat), nrow = 1223, ncol = 60660, byrow = FALSE,
                        dimnames = dimnames(bulk_rna_mat))
saveRDS(bulk_rna_mat2, "bulk_matrix.rds")

# getting cell type data
sc_meta <- read.csv("metadata_clean.csv", header = T)
cell.type.labels <- sc_meta[sc_meta$cell_ids %in% row.names(scrnaseq_mat), "celltype_major"]
# getting cell state labels for only cells in cell_types
cell.state.labels <- sc_meta[sc_meta$cell_ids %in% row.names(scrnaseq_mat), "cellType"]

## QC of cell states and cell labels
plot.cor.phi (input=scrnaseq_mat,
              input.labels=cell.state.labels,
              title="cell state correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.cs", 
              cexRow=0.2, cexCol=0.2,
              margins=c(2,2))
plot.cor.phi (input=scrnaseq_mat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              cexRow=0.5, cexCol=0.5,
)

# detecting outliers
sc.stat <- plot.scRNA.outlier(
  input=scrnaseq_mat,  
  cell.type.labels=cell.type.labels,
  species="hs", 
  return.raw=TRUE 
  )
write.csv(sc.stat, "sc_outliers.csv")
# bk.stat is giving error
bk.stat <- plot.bulk.outlier(
  bulk.input=bulk_rna_mat2, 
  sc.input=scrnaseq_mat,  
  cell.type.labels=cell.type.labels,
  species="hs", 
  return.raw=TRUE
  )

# column names from stats are the gene groups with outliers
# ribosomal protein, mt ribosomal protein, ribosomal pseudogenes and
# chrM were detected as outliers in both datasets
sc.dat.filtered <- cleanup.genes (input=scrnaseq_mat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=colnames(sc.stat)[3:11],
                                  exp.cells=5)
#saveRDS(sc.dat.filtered, "sc_no_outliers.rds")
# correlation between bulk and filtered scrnaseq
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                bulk.input = bulk_rna_mat2)
# using protein coding genes
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")
saveRDS(sc.dat.filtered.pc, "proteinCoding.rds")
# 14328 genes selected
# rm("sc.dat.filtered")
# differential expression to select marker genes
diff.exp.stat <- get.exp.stat(sc.dat=scrnaseq_mat[,colSums(scrnaseq_mat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=16 #number of threads
)
saveRDS(diff.exp.stat, "diff_expression.rds")
# generating count matrix for marker genes alone
sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)
saveRDS(sc.dat.filtered.pc.sig, "filtered_withmarker.rds")
## to create prism object on filtered signficant genes data generated above.

# creating prism object
wu_etal_prism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bulk_rna_mat2,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL, # we dont have cell labelled as "tumor", all cells treated in similar way
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
saveRDS(wu_etal_prism, "brca_prism.rds")
bp.res <- run.prism(prism = wu_etal_prism, n.cores=16) # too low cores, ~30 min runtime

# extract cell type fraction
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
saveRDS(theta, "theta_data.rds")
theta.cv <- bp.res@posterior.theta_f@theta.cv
save(bp.res, file="wu_etal_deconvoluted_results.rds")

theta.df <- as.data.frame(theta)
theta.cv.df <- as.data.frame(theta.cv)
plot(x=theta.cv.df$Endothelial, theta.df$Endothelial) # negative correlation
write.csv(theta.cv.df, "cell_deconvolution.csv")
write.csv(theta.df, "cell_deconvolution_fullDF.csv")


