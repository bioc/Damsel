#' Plot correlation heatmap
#'
#' `plotCorrHeatmap` plots the correlation of all available BAM files Dam and Fusion, to visualise the similarity between files.
#' * uses the non-parametric "spearman's" correlation.
#'
#' The correlation between Dam_1 and Fusion_1 can be expected to reach ~ 0.7, whereas the correlation between Dam_1 & Dam_3 or Fusion_1 & Fusion_2 would be expected to be closer to ~0.9
#'
#' @param df A data.frame of GATC region counts as outputted from [countBamInGatc()].
#' @param method The correlation method used. If not specified, will use default of non-parametric spearman's.
#' * Non-parametric methods are recommended as data does not reliably meet the requirements for parametric analysis.
#'
#' @return A heatmap style plot of the samples, coloured by correlation value.
#' * Colour spectrum is determined from the minimum correlation as the lowest correlation, the median correlation as the midpoint colour, and 1 as the top colour.
#' @examples
#' counts.df <- random_counts()
#' plotCorrHeatmap(counts.df, method = "spearman")
#' @export
plotCorrHeatmap <- function(df, method="spearman") {
    if (!is.data.frame(df)) {
        stop("data.frame of counts is required")
    }
    if (missing(method)) {
        message("default spearman's method is used")
    }

    corr_res <- stats::cor(df[, grepl("bam", colnames(df), ignore.case = TRUE)], method = method)
    ComplexHeatmap::Heatmap(corr_res, name = "Colour Legend")
}
