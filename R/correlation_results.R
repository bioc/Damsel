#' Plot correlation heatmap
#'
#' `plotCorrHeatmap` plots the correlation of all available BAM files Dam and Fusion, to visualise the similarity between files.
#' * uses the non-parametric "spearman's" correlation.
#'
#' The correlation between Dam_1 and Fusion_1 can be expected to reach ~ 0.7, whereas the correlation between Dam_1 & Dam_3 or Fusion_1 & Fusion_2 would be expected to be closer to ~0.9
#'
#' @param df A data.frame of GATC region counts as outputted from [process_bams()].
#' @param method The correlation method used. If not specified, will use default of non-parametric spearman's.
#' * Non-parametric methods are recommended as data does not reliably meet the requirements for parametric analysis.
#'
#' @return A `ggplot2` object. A heatmap style plot of the samples, coloured by correlation value.
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
    median_corr <- round(stats::median(corr_res), 1)
    min_corr <- floor(min(corr_res) * 10) / 10
    corr_res <- round(corr_res, 2)
    # Use correlation between variables as distance and reorder
    dd <- stats::as.dist((1 - corr_res) / 2)
    hc <- stats::hclust(dd)
    corr_res <- corr_res[hc$order, hc$order]
    # upper triangle
    corr_res[lower.tri(corr_res)] <- NA
    # Melt the correlation matrix
    corr_res <- reshape2::melt(corr_res, na.rm = TRUE)
    # plot heatmap
    ggplot2::ggplot(corr_res, ggplot2::aes(.data$Var2, .data$Var1, fill = .data$value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_gradient2(
            low = "blue", high = "red",
            midpoint = median_corr, limit = c(min_corr, 1), space = "Lab",
            name = paste(stringr::str_to_title(method), "\nCorrelation", sep = " ")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
        ggplot2::coord_fixed()
}
