#' Create DGE object for differential testing
#'
#' `edgeR_set_up()` sets up the edgeR analysis for visualisation of the samples [edgeR_plot_mds()], and then for identifying differentially methylated regions [edgeR_results()].
#'
#'
#' @param df data.frame generated from [process_bams]. Ensure that the samples are ordered by (Dam_1.bam, Fusion_1.bam, Dam_2.bam, Fusion_2.bam, ...)
#' @param lib.size Library size for each sample is calculated as the sum across all rows for that sample unless otherwise specified
#' @param keep_a Minimum cpm of the counts - default is 0.5
#' @param keep_b Minimum number of samples to meet the criteria of keep_a in order to retain the region in the downstream analysis. Default is 3 (assuming 6 samples)
#'
#' @return Refer to [edgeR::?`DGEListClass`] for details
#' @export
#'
#' @examples
#' counts.df <- random_counts()
#'
#' edgeR_set_up(counts.df)
#dmSetUp
edgeR_set_up <- function(df, lib.size=NULL, keep_a=0.5, keep_b=3) {
  if(!is.data.frame(df)) {
    stop("Must have data.frame of counts")
  }
  if(!is.numeric(keep_a) | length(keep_a) > 1) {
    stop("keep_a must be 1 value, recommend using default value")
  }
  if(!is.numeric(keep_b) | length(keep_b) > 1) {
    stop("keep_b must be 1 value, recommend using default value")
  }

  matrix <- as.matrix(df[,grepl("bam", colnames(df))]) # can I be sure they would have "bam" in it?
  rownames(matrix) <- df$Position

  n_samples <- seq(1:(ncol(matrix)/2))

  group <- rep(c("Dam", "Fusion"), times = length(n_samples)) #specify that this order is required

  dge <- edgeR::DGEList(matrix, lib.size = lib.size, group = group, gene = df[,2:5])

  keep <- rowSums(edgeR::cpm(dge) >= keep_a) >= keep_b #potentially will mess with
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  dge <- edgeR::calcNormFactors(dge)

  design_df <- seq(1:ncol(matrix)) %>%
                 data.frame() %>%
                 stats::setNames("group")
  zero_vec <- rep(0, times = length(n_samples))
  #need to replace the 1:length with seq along or something
  for(i in 1:length(n_samples)) {
    design <- replace(zero_vec, i, 1)
    design <- rep(design, each = 2)
    design_df[,ncol(design_df) + 1] <- design

  }
  design <- stats::model.matrix(~., data = design_df[, 1:ncol(design_df) - 1])

  dge <- edgeR::estimateDisp(dge, robust = T, design = design)

  dge
}

#' Plot differences between samples: Multidimensional scaling
#'
#' `edgeR_plot_mds` visualises the difference between samples.
#' * expect control (Dam-only) samples to cluster together and for Fusion samples to cluster together
#'
#' @param dge as outputted from [edgeR_set_up()]
#'
#' @return refer to [edgeR::plotMDS()] for details
#' @export
#'
#' @examples
#' counts.df <- random_counts()
#' dge <- edgeR_set_up(counts.df)
#'
#' edgeR_plot_mds(dge)
#dmPlotMDS
edgeR_plot_mds <- function(dge) {
  group <- dge$samples$group %>% as.character()
  limma::plotMDS(dge, col = as.numeric(factor(group)))
}


#' Differential testing
#'
#' `edgeR_results` calculates the differential methylation results, identifying which GATC regions have been enriched in the Fusion samples relative to the controls.
#' Refer to the following pages for further details:
#' * [edgeR::glmQLFit()]
#' * [edgeR::glmGQLFTest()]
#' * [edgeR::decideTestsDGE()]
#'
#' @param dge as outputted from [edgeR_set_up()]
#' @param p.value p value threshold for minimum significance. Default is 0.05
#' @param lfc minimum log fold change for significant results. Default is 1
#' @param regions data.frame of GATC regions. If not provided, default used is `regions_gatc_drosophila_dm6`
#'
#' @return data.frame of differential methylation results.
#' Columns are as follows;
#' * rownames(Region position),
#' * logFC (log fold change),
#' * logCPM (log counts per million),
#' * F (F statistic used to identify significance),
#' * PValue,
#' * adjustedP (P Value with multiple test correction),
#' * de (result: -1,0,1)
#' @export
#'
#' @examples
#' counts.df <- random_counts()
#' dge <- edgeR_set_up(counts.df)
#'
#' edgeR_results <- edgeR_results(dge, p.value = 0.05, lfc = 1)
#' head(edgeR_results)
#dmResults
#also need to update this fn - adjusted p val
edgeR_results <- function(dge, p.value=0.05, lfc=1, regions=regions_gatc_drosophila_dm6) {
  if(!is.numeric(p.value) | length(p.value) > 1) {
    stop("p.value must be 1 number, recommend using default value")
  }
  if(!is.numeric(lfc) | length(lfc) > 1) {
    stop("lfc must be 1 number, recommend using default value")
  }
  if(missing(regions)) {
    message("regions missing, default `regions_gatc_drosophila_dm6` used instead")
  }
  group <- dge$samples$group %>% as.character()
  design <- dge$design
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)
  #summary(de <- edgeR::decideTestsDGE(qlf, p.value = p.value, lfc = lfc))
  lrt_table <- qlf$table
  lrt_table <- lrt_table %>% dplyr::mutate(adjust.p = stats::p.adjust(PValue, method = "BH"),
                                    de = dplyr::case_when(logFC < -lfc & adjust.p < p.value ~ -1,
                                                   abs(logFC) < lfc ~ 0,
                                                   logFC > lfc & adjust.p < p.value ~ 1, TRUE ~ 0))
  lrt_table <- add_de(lrt_table, regions)
  lrt_table
}


#' Plot differential testing results
#'
#' `edgeR_results_plot` provides an MA style plot for differential methylation results, allowing for a visualisation of the logFC, P values, and spread of -1,0,1 results.
#' * for further details, see [edgeR::plotSmear()]
#'
#' @param dge as outputted from [edgeR_set_up()]
#' @param p.value p value threshold for minimum significance. Default is 0.05
#' @param lfc minimum log fold change for significant results. Default is 1
#'
#' @return MA style scatter plot with average logCPM on x-axis, average logFC on y-axis, with dots coloured by significance
#' @export
#'
#' @examples
#' counts.df <- random_counts()
#' dge <- edgeR_set_up(counts.df)
#'
#' edgeR_results_plot(dge, p.value = 0.05, lfc = 1)
#dmPlotResults
edgeR_results_plot <- function(dge, p.value=0.05, lfc=1) {
  if(!is.numeric(p.value) | length(p.value) > 1) {
    stop("p.value must be 1 number, recommend using default value")
  }
  if(!is.numeric(lfc) | length(lfc) > 1) {
    stop("lfc must be 1 number, recommend using default value")
  }
  group <- dge$samples$group %>% as.character()
  design <- dge$design
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)
  #summary(de <- edgeR::decideTestsDGE(qlf, p.value = p.value, lfc = lfc))
  lrt_table <- qlf$table
  lrt_table <- lrt_table %>% dplyr::mutate(adjust.p = stats::p.adjust(PValue, method = "BH"),
                                           de = dplyr::case_when(logFC < -lfc & adjust.p < p.value ~ -1,
                                                                 abs(logFC) < lfc ~ 0,
                                                                 logFC > lfc & adjust.p < p.value ~ 1, TRUE ~ 0))
  detags <- rownames(dge)[as.logical(lrt_table$de)]
  edgeR::plotSmear(qlf, de.tags=detags, ylab = "logFC - Fusion/Dam")
}

