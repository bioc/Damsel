#' Plot distribution of counts
#' `plotCountsDistribution` plots the distribution of counts enabling the comparison of different samples.
#'  Can highlight the different library sizes of the samples.
#' @param counts.df A counts data.frame as outputted from `countBamInGatc`
#' @param constant A numerical offset to avoid log2(0) = -Inf. Default is 1
#'
#' @return A ggplot2 density plot
#' @export
#'
#' @examples
#' set.seed(123)
#' counts.df <- random_counts()
#'
#' plotCountsDistribution(counts.df, 1)
plotCountsDistribution <- function(counts.df, constant=1) {
    df <- counts.df
    plot_df <- df[,grepl("bam", colnames(df), ignore.case = TRUE)]
    if(ncol(plot_df) == 0) {
      stop("No counts identified in data.frame. Please use counts.df with .bam/.BAM extension")
    }
    plot_df$Position <- df$Position
    plot_df <- plot_df %>% reshape2::melt(id = "Position")

    plot <- plot_df %>%
        ggplot2::ggplot(ggplot2::aes(x = log2(.data$value + constant),
            colour = .data$variable)) +
        ggplot2::geom_density() +
        ggplot2::labs(x = paste0("Log2(Count + ", constant, ")"), y = "Density")
    plot +
        ggplot2::theme_classic()
}



#' Plotting the
#'
#' @param counts.df A counts data.frame as outputted from `countBamInGatc`
#' @param dm_results.df A data.frame of differential testing results as outputted from `testDmRegions()`.
#' @param peaks.df A data.frame of peaks as outputted from `identifyPeaks()`.
#' @param position If the bar plots should be stacked (showing the count), or fill (showing proportion). Default is "stack".
#'
#' @return A ggplot2 bar plot
#' @export
#'
#' @examples
#' set.seed(123)
#' counts.df <- random_counts()[1:2,]
#' dm_results <- random_edgeR_results()[1:2,]
#' peaks <- identifyPeaks(dm_results)
#'
#' # stacked plot
#' plotCountsInPeaks(counts.df, dm_results, peaks, position = "stack")
#' # filled plot
#' plotCountsInPeaks(counts.df, dm_results, peaks, position = "fill")
plotCountsInPeaks <- function(counts.df, dm_results.df, peaks.df,
    position=c("stack", "fill")) {
    if (is.null(counts.df) | is.null(dm_results.df) |
        is.null(peaks.df)) {
        stop("All data.frames must be inputted (counts.df, dm_results.df, peaks.df)")
    }
    dm_results <- dplyr::filter(dm_results.df, Position %in% counts.df$Position)
    counts <- counts.df %>%
        dplyr::mutate(tested = dm_results$meth_status,
            tested = ifelse(.data$tested == "Not_included", FALSE, TRUE))
    overlap <- plyranges::find_overlaps(plyranges::as_granges(counts),
        plyranges::as_granges(peaks.df)) %>%
        data.frame()
    df <- counts[,grepl("bam", colnames(counts), ignore.case = TRUE)]
    if(nrow(overlap) == 0) {
        df <- dplyr::mutate(df, in_peak = FALSE)
    } else {
        df <- df %>%
            dplyr::mutate(in_peak = ifelse(counts$Position %in% overlap$Position, TRUE, FALSE))
    }
    df <- df %>%
        dplyr::mutate(in_peak = ifelse(counts$tested == FALSE, "Not_tested", .data$in_peak)) %>%
        reshape2::melt(id = "in_peak")
    if(length(position) > 1) {
        position <- "stack"
    }
    plot <- df %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$variable, y = .data$value,
            fill = .data$in_peak)) +
        ggplot2::geom_col(position = position) +
        ggplot2::labs(x = "Sample", y = "Count")
    plot +
        ggplot2::theme_bw()
}

