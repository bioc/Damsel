#' Create DGE object for differential testing
#'
#' `edgeR_set_up()` sets up the edgeR analysis for visualisation of the samples [edgeR_plot_mds()], and then for identifying differentially methylated regions [edgeR_results()].
#'
#' @param counts.df A data.frame generated from [process_bams]. Ensure that the samples are ordered by (Dam_1.bam, Fusion_1.bam, Dam_2.bam, Fusion_2.bam, ...)
#' @param max.width Remove large regions, default is width of 10,000. We recommend this value as the Dam can methylate GATC sites up to 5kb away from the binding site, generating a total width of 10 kb.
#' @param lib.size Library size for each sample is calculated as the sum across all rows for that sample unless otherwise specified
#' @param keep_a Filtering parameter, minimum counts per million (cpm) of each sample. Recommend leaving at default of 0.5
#' @param keep_b Filtering parameter, minimum number of samples to meet the criteria of keep_a in order to retain the region in the downstream analysis. Default is 3 (assuming 6 samples)
#'
#' @return An object of class `DGEList`. Refer to [edgeR::?`DGEListClass`] for details
#' @export
#'
#' @examples
#' counts.df <- random_counts()
#'
#' edgeR_set_up(counts.df)
#dmSetUp
edgeR_set_up <- function(counts.df, max.width=10000, lib.size=NULL, keep_a=0.5, keep_b=3) {
  if(!is.data.frame(counts.df)) {
    stop("Must have data.frame of counts")
  }
  if(!is.numeric(keep_a) | length(keep_a) > 1) {
    stop("keep_a must be 1 value, recommend using default value")
  }
  if(!is.numeric(keep_b) | length(keep_b) > 1) {
    stop("keep_b must be 1 value, recommend using default value")
  }

  counts.df <- counts.df %>% dplyr::filter(.data$width <= max.width)
  matrix <- as.matrix(counts.df[,grepl("bam", colnames(counts.df))]) # can I be sure they would have "bam" in it?
  rownames(matrix) <- counts.df$Position

  n_samples <- seq_len(ncol(matrix)/2)

  group <- rep(c("Dam", "Fusion"), times = length(n_samples)) #specify that this order is required

  dge <- edgeR::DGEList(matrix, lib.size = lib.size, group = group, gene = counts.df[,2:5])

  keep <- rowSums(edgeR::cpm(dge) >= keep_a) >= keep_b #potentially will mess with
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  dge <- edgeR::calcNormFactors(dge)

  design_df <- seq_len(ncol(matrix)) %>%
                 data.frame() %>%
                 stats::setNames("group")
  zero_vec <- rep(0, times = length(n_samples))
  #need to replace the 1:length with seq along or something
  for(i in n_samples) {
    design <- replace(zero_vec, i, 1)
    design <- rep(design, each = 2)
    design_df[,ncol(design_df) + 1] <- design

  }
  design <- stats::model.matrix(~., data = design_df[, seq_len(ncol(design_df)) - 1])

  dge <- edgeR::estimateDisp(dge, robust = TRUE, design = design)

  dge
}

#' Plot differences between samples: Multidimensional scaling
#'
#' `edgeR_plot_mds` visualises the difference between samples.
#' * expect control (Dam-only) samples to cluster together and for Fusion samples to cluster together
#'
#' @param dge A DGEList object as outputted from [edgeR_set_up()]
#'
#' @return An object of class `MDS`. Refer to [edgeR::plotMDS()] for details
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
#' * [edgeR::glmQLFTest()]
#' * [edgeR::decideTestsDGE()]
#'
#' @param dge A DGEList object as outputted from [edgeR_set_up()]
#' @param p.value A number between 0 and 1 providing the required false discovery rate (FDR). Default is 0.01
#' @param lfc A number giving the minimum absolute log2-fold-change for significant results. Default is 1
#' @param regions A data.frame of GATC regions. If not provided, default used is `regions_gatc_drosophila_dm6`
#'
#' @return A `data.frame` of differential methylation results. Columns are: Position (chromosome-start), seqnames, start, end, width, strand, number (region number), dm (edgeR result: -1,0,1,NA), logFC, adjust.p, meth_status (Downreg, No_signal, Upreg, Not_included)
#' @export
#'
#' @examples
#' set.seed(123)
#' example_regions <- random_regions()
#' counts.df <- random_counts()
#' dge <- edgeR_set_up(counts.df)
#'
#' dm_results <- edgeR_results(dge, p.value = 0.01, lfc = 1, regions = example_regions)
#' head(dm_results)
#dmResults
#also need to update this fn - adjusted p val
edgeR_results <- function(dge, p.value=0.01, lfc=1, regions) {
  if(!is.numeric(p.value) | length(p.value) > 1) {
    stop("p.value must be 1 number, recommend using default value")
  }
  if(!is.numeric(lfc) | length(lfc) > 1) {
    stop("lfc must be 1 number, recommend using default value")
  }
  if(missing(regions)) {
    message("GATC region data.frame required")
  }
  group <- dge$samples$group %>% as.character()
  design <- dge$design
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)
  lrt_table <- qlf$table
  lrt_table <- lrt_table %>% dplyr::mutate(adjust.p = stats::p.adjust(.data$PValue, method = "BH"),
                                    dm = dplyr::case_when(.data$logFC < lfc & .data$adjust.p < p.value ~ -1,
                                                   abs(.data$logFC) < lfc ~ 0,
                                                   .data$logFC > lfc & .data$adjust.p < p.value ~ 1, TRUE ~ 0))
  lrt_table <- add_de(lrt_table, regions)
  lrt_table
}


#' Plot differential testing results
#'
#' `edgeR_results_plot` provides an MA style plot for differential methylation results, allowing for a visualisation of the logFC, P values, and spread of -1,0,1 results.
#' * for further details, see [edgeR::plotSmear()]
#'
#' @param dge A DGEList object as outputted from [edgeR_set_up()]
#' @param p.value A number between 0 and 1 providing the required false discovery rate (FDR). Default is 0.01
#' @param lfc A number giving the minimum absolute log2-fold-change for significant results. Default is 1
#'
#' @return An object of class `MA` style scatter plot with average logCPM on x-axis, average logFC on y-axis, with dots coloured by significance
#' @export
#'
#' @examples
#' counts.df <- random_counts()
#' dge <- edgeR_set_up(counts.df)
#'
#' edgeR_results_plot(dge, p.value = 0.01, lfc = 1)
#dmPlotResults
edgeR_results_plot <- function(dge, p.value=0.01, lfc=1) {
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
  lrt_table <- qlf$table
  lrt_table <- lrt_table %>% dplyr::mutate(adjust.p = stats::p.adjust(.data$PValue, method = "BH"),
                                           dm = dplyr::case_when(.data$logFC < -lfc & .data$adjust.p < p.value ~ -1,
                                                                 abs(.data$logFC) < lfc ~ 0,
                                                                 .data$logFC > lfc & .data$adjust.p < p.value ~ 1, TRUE ~ 0))
  detags <- rownames(dge)[as.logical(lrt_table$dm)]
  edgeR::plotSmear(qlf, de.tags=detags, ylab = "logFC - Fusion/Dam")
}


#' Add dm_results to GATC regions
#'
#' Used within aggregate_peaks to add in the regions that were excluded from edgeR analysis for low counts
#'
#' @param dm_results A data.frame as outputted from [edgeR_results()]
#' @param regions A data.frame of GATC regions, default is regions_gatc_drosophila_dm6
#'
#' @return A `data.frame` of regions with added information about the dm results;
#' * dm - 1,0,-1,NA ;
#' * logFC: 0 if dm is NA ;
#' * adjust.p: 1 if dm is NA :
#' * meth_status: Upreg, No_sig, Downreg, Not_included
add_de <- function(dm_results, regions) {
  if(!is.data.frame(dm_results)) {
    stop("Must have data frame of differential testing results from `edgeR_results")
  }
  if(!is.data.frame(regions)) {
    stop("Regions must be a data.frame")
  }
  results <- dm_results
  df <- regions %>%
    dplyr::mutate(seqnames = paste0("chr", .data$seqnames), number = seq_len(nrow(.)))
  df$dm <- results[match(df$Position, row.names(results)), "dm"]
  df$logFC <- results[match(df$Position, row.names(results)), "logFC"]
  df$PValue <- results[match(df$Position, row.names(results)), "PValue"]
  df$adjust.p <- results[match(df$Position, row.names(results)), "adjust.p"]
  df <- df %>%
    dplyr::mutate(meth_status = dplyr::case_when(is.na(.data$dm) ~ "Not_included",
                                                 .data$dm == 1 ~ "Upreg",
                                                 .data$dm == -1 ~ "Downreg",
                                                 TRUE ~ "No_sig"))
  df <- df %>%
    dplyr::mutate(logFC = dplyr::coalesce(.data$logFC, 0),
                  PValue = dplyr::coalesce(.data$PValue, 1),
                  adjust.p = dplyr::coalesce(.data$adjust.p, 1))
  df
}
