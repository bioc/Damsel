#' Plotting peaks
#'
#' `geom_peak.new` is a ggplot layer that visualises the positions of peaks across a given region.
#' * cannot be plotted by itself, must be added to an existing ggplot object - see examples.
#'
#'
#' @param peaks.df A data.frame of peaks as outputted from `aggregate_peaks()`
#' @param peak.label Specify whether peak_id labels should be added to the plot. Default is FALSE
#' @param peak.color Specify colour of peak. Default is black
#' @param peak.size Specify size of rectangle. Default is 5
#' @param plot.space Specify gap to next plot. Recommend leaving to the default: 0.1
#' @param plot.height Specify overall height of plot. Recommend leaving to the default: 0.05
#'
#' @return A `ggplot_add` object.
#' @export
#'
#' @examples
#' set.seed(123)
#' counts.df <- random_counts()
#' dm_results <- random_edgeR_results()
#' peaks <- aggregate_peaks(dm_results)
#' plot_counts_all_bams(counts.df,
#'     seqnames = "chr2L",
#'     start_region = 1,
#'     end_region = 40000,
#'     log2_scale = FALSE
#' ) +
#'     geom_peak.new(peaks)
#'
#' plot_counts_all_bams(counts.df,
#'     seqnames = "chr2L",
#'     start_region = 1,
#'     end_region = 40000,
#'     log2_scale = FALSE
#' ) +
#'     geom_peak.new(peaks, peak.label = TRUE)
#' # The plots can be layered -------------------------------------------------
geom_peak.new <- function(peaks.df = NULL, peak.label = FALSE, peak.color = "black", peak.size = 5,
                          plot.space = 0.1, plot.height = 0.05) {
    structure(
        list(
            peaks.df = peaks.df, peak.label = peak.label, peak.color = peak.color, peak.size = peak.size,
            plot.space = plot.space, plot.height = plot.height
        ),
        class = "peak.new"
    )
}


#' @export
ggplot_add.peak.new <- function(object, plot, object_name) {
    if (!is.data.frame(object$peaks.df)) {
        stop("data.frame of peaks is required")
    }
    plot2 <- plot
    while ("patchwork" %in% class(plot2)) {
        plot2 <- plot2[[1]]
    }
    plot.data <- plot2$labels$title
    plot.data <- stringr::str_split_1(plot.data, ":")

    plot.chr <- plot.data[1]
    plot.data <- stringr::str_split_1(plot.data[2], "-")
    plot.region.start <- plot.data[1] %>% as.numeric()
    plot.region.end <- plot.data[2] %>% as.numeric()

    # get parameters
    peaks.df <- object$peaks.df
    peak.label <- object$peak.label
    peak.color <- object$peak.color
    peak.size <- object$peak.size
    plot.space <- object$plot.space
    plot.height <- object$plot.height

    valid.bed <- GetRegion_hack(df = peaks.df, columns = c("seqnames", "start", "end", "peak_id"), chr = plot.chr, start = plot.region.start, end = plot.region.end)

    peak.plot <- plot_peak(valid.bed = valid.bed, plot.size = peak.size, plot.color = peak.color, peak.label = peak.label)

    peak.plot <- peak.plot + ggplot2::labs(y = "Peak") +
        theme_peak_hack(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end))
    # assemble plot
    patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
        peak.plot,
        ncol = 1, heights = c(1, plot.height)
    )
}
