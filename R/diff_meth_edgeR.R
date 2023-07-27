
#' Differential methylation analysis
#'
#' `edgeR_set_up()` sets up the edgeR analysis for visualisation of the samples [edgeR_plot_mds()], and then for identifying differentially methylated regions [edgeR_results()].
#'
#'
#' @param df Data frame generated from [process_bams]. Ensure that the samples are ordered by (Dam_1, Fusion_1, Dam_2, Fusion_2, ...)
#' @param lib.size Library size for each sample is calculated as the sum across all rows for that sample unless otherwise specified
#' @param keep_a Minimum cpm of the counts - default is 0.5
#' @param keep_b Minimum number of samples to meet the criteria of keep_a in order to retain the region in the downstream analysis. Default is 3 (assuming 6 samples)
#'
#' @return Refer to [edgeR::?`DGEListClass`] for details
#' @export
#'
#' @examples
edgeR_set_up <- function(df, lib.size=NULL, keep_a=0.5, keep_b=3) {
  matrix <- as.matrix(df[,grepl("bam", colnames(df))]) # can I be sure they would have "bam" in it?
  rownames(matrix) <- df$Position

  n_samples <- seq(1:(ncol(test_peaks_bams)/2))

  group <- rep(c("Dam", "Fusion"), times = length(n_samples)) #specify that this order is required

  dge <- edgeR::DGEList(matrix, lib.size = lib.size, group = group, gene = df[,2:5])

  keep <- rowSums(cpm(dge) >= keep_a) >= keep_b #potentially will mess with
  dge <- dge[keep, , keep.lib.sizes=FALSE]

  dge <- edgeR::calcNormFactors(dge)

  design_df <- seq(1:10) %>% as.data.frame() %>% setNames("group")
  zero_vec <- rep(0, times = length(n_samples))
  for(i in 1:length(n_samples)) {
    design <- replace(zero_vec, i, 1)
    design <- rep(design, each = 2)
    design_df[,ncol(design_df) + 1] <- design

  }
  design <- model.matrix(~., data = design_df[,1:ncol(design_df)-1])

  dge <- edgeR::estimateDisp(dge, robust = T, design = design)

  dge
}

#' Multidimensional scaling plot for samples
#'
#' this is a test for having 2 in 1 doc
#' @param dge as outputted from [edgeR_set_up()]
#'
#' @return mds plot of samples -
#' @export
#'
#' @examples
edgeR_plot_mds <- function(dge) {
  group <- dge$samples$group %>% as.character()
  edgeR::plotMDS(dge, col=as.numeric(factor(group)))
}

edgeR_results <- function(dge, p.value=0.05, lfc=1) {
  group <- dge$samples$group %>% as.character()
  design <- dge$design
  fit = edgeR::glmQLFit(dge, design = design)
  qlf = edgeR::glmQLFTest(fit, coef = 2)
  summary(de <- edgeR::decideTestsDGE(qlf, p.value = p.value, lfc = lfc))
  lrt_table <- qlf$table
  lrt_table$de <- data.frame(de)$groupFusion
  write.table(lrt_table, file='../output/lrt_sd.txt', quote=F)
  write.table(keep, file='../output/keep', quote=F, col.names = FALSE)
  lrt_table
}

edgeR_results_plot <- function(dge, results) {
  detags <- rownames(dge)[as.logical(results$de)]
  edgeR::plotSmear(qlf, de.tags=detags, ylab = "logFC - Scalloped/Dam")
}

