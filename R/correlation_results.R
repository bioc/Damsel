#' Plotting correlation of all samples
#'
#' `corr_heatmap` plots the correlation of all available BAM files Dam and Fusion, to visualise the similarity between files.
#' * uses the non-parametric "spearman's" correlation.
#'
#' The correlation between Dam_1 and Fusion_1 can be expected to reach ~ 0.7, whereas the correlation between Dam_1 & Dam_3 or Fusion_1 & Fusion_2 would be expected to be closer to ~0.9
#'
#' @param df data frame as outputted from [process_bams()]
#' @param method correlation method. Default is the non-parametric spearman's. Non-parametric methods are recommended as data does not reliably meet the requirements for parametric analysis.
#'
#' @return A heatmap style plot of the samples, coloured by correlation value. Colour spectrum is determined from the minimum correlation as the lowest correlation, the median correlation as the midpoint colour, and 1 as the top colour.
#' @examples
#'
#' @export
#corrHeatmap
corr_heatmap <- function(df, method="spearman") {
  if(missing(df) | !is.data.frame(df)) {
    stop("data.frame of counts is required")
  }
  if(missing(method)) {
    message("default spearman's method is used")
  }

  corr_res <- cor(df[, grepl("bam", colnames(df))], method = method)
  median_corr <- round(median(corr_res), 1)
  min_corr <- floor(min(corr_res) * 10) / 10
  corr_res <- round(corr_res, 2)
  # Use correlation between variables as distance and reorder
  dd <- as.dist((1 - corr_res) / 2)
  hc <- hclust(dd)
  corr_res <- corr_res[hc$order, hc$order]
  # upper triangle
  corr_res[lower.tri(corr_res)] <- NA
  # Melt the correlation matrix
  corr_res <- reshape2::melt(corr_res, na.rm = TRUE)
  #plot heatmap
  ggplot2::ggplot(corr_res, ggplot2::aes(Var2, Var1, fill = value)) +
    ggplot2::geom_tile(color = "white")+
    ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                                 midpoint = median_corr, limit = c(min_corr,1), space = "Lab",
                                 name = paste(stringr::str_to_title(method), "\nCorrelation", sep = " ")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
    ggplot2::coord_fixed()
}

#' Plotting correlation of specific samples
#'
#' `corr_scatter` plots a scatterplot with correlation results from two selected samples, in order to visualise their similarity/difference.
#' * uses the non-parametric "spearman's" correlation.
#' * for more information, see: [ggpubr::ggscatter()]
#'
#'The correlation between Dam_1 and Fusion_1 can be expected to reach ~ 0.7, whereas the correlation between Dam_1 & Dam_3 or Fusion_1 & Fusion_2 would be expected to be closer to ~0.9
#'
#' @param df data frame as outputted from [process_bams()]
#' @param sample_1 name of the BAM file within df
#' @param sample_2 name of the BAM file within df
#' @param method correlation method. Default is the non-parametric spearman's. Non-parametric methods are recommended as data does not reliably meet the requirements for parametric analysis.
#'
#'
#' @return A scatterplot of the two selected samples, showing their counts at each region, overlaid with the correlation results.
#' @export
#'
#' @examples
#corrScatter
corr_scatter <- function(df, sample_1, sample_2, method="spearman") {
  if(missing(df) | !is.data.frame(df)) {
    stop("data.frame of counts is required")
  }
  if(missing(sample_1) | !is.character(sample_1)) {
    stop("sample_1 must be a character vector")
  }
  if(missing(sample_2) | !is.character(sample_2)) {
    stop("sample_2 must be a character vector")
  }
  if(missing(method)) {
    message("default spearman's method is used")
  }

  ggpubr::ggscatter(df, x = sample_1, y = sample_2,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = TRUE, cor.method = method)
}
