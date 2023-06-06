
edgeR_set_up <- function(df, lib.size=NULL, keep_a=0.5, keep_b=3) {
  matrix <- as.matrix(df[,7:ncol(df)])
  rownames(matrix) <- df$Position

  group <- c("Dam", "Fusion", "Dam", "Fusion", "Dam", "Fusion") #how can I ensure this order?

  dge <- edgeR::DGEList(matrix, lib.size = lib.size, group = group, gene = df[,2:5])

  keep <- rowSums(cpm(dge) >= keep_a) >= keep_b #potentially will mess with
  dge <- dge[keep, , keep.lib.sizes=FALSE]

  dge <- edgeR::calcNormFactors(dge)

  design <- model.matrix(~group + c(1,1,0,0,0,0) + c(0,0,1,1,0,0)) #depends on N of samples

  dge <- edgeR::estimateDisp(dge, robust = T, design = design)

  dge
}

edgeR_plot_mds <- function(dge, group) {
  edgeR::plotMDS(dge, col=as.numeric(factor(group))) #work on this
}

edgeR_results <- function(dge, design, p.value=0.05, lfc=1) {
  fit = edgeR::glmQLFit(dge, design = design)
  qlf = edgeR::glmQLFTest(fit, coef = 2)
  summary(de <- edgeR::decideTestsDGE(qlf, p.value = p.value, lfc = lfc))
  lrt_table <- qlf$table
  lrt_table$de <- as.data.frame(de)$groupFusion
  write.table(lrt_table, file='../output/lrt_sd.txt', quote=F)
  write.table(keep, file='../output/keep', quote=F, col.names = FALSE)
  lrt_table
}

edgeR_results_plot <- function(dge, de) {
  detags <- rownames(dge)[as.logical(de)]
  edgeR::plotSmear(qlf, de.tags=detags, ylab = "logFC - Scalloped/Dam")
}

