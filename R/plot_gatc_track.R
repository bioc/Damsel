#' Plot for a GATC track
#'
#' `geom_gatc` is a ggplot layer that visualises the positions of GATC sites across a given region.
#' * cannot be plotted by itself, must be added to an existing ggplot object - see examples.
#'
#' @param gatc_sites.df A data.frame of positions of GATC sites - can be made from `gatc_track()$sites`
#' @param gatc.color Specify colour of lines. Default is red
#' @param gatc.size Specify size of the line. Default is 5
#' @param plot.space Specify gap to next plot. Recommend leaving to the default: 0.2
#' @param plot.height Specify overall height of the plot. Recommend leaving to the default: 0.05
#'
#' @return A `ggplot_add` object.
#' @export
#' @references ggcoverage
#' @seealso [geom_peak()] [geom_dm()] [geom_genes()] [geom_counts()] [plot_wrap()]
#' @examples
#' set.seed(123)
#' example_regions <- random_regions()
#' counts.df <- random_counts()
#' gatc_sites <- dplyr::mutate(example_regions,
#'     seqnames = paste0("chr", seqnames),
#'     start = start - 3, end = start + 4
#' )
#'
#' plot_counts_all_bams(counts.df,
#'     seqnames = "chr2L",
#'     start_region = 1,
#'     end_region = 40000,
#'     log2_scale = FALSE
#' ) +
#'     geom_gatc(gatc_sites)
#' # The plots can be layered -------------------------------------------------
geom_gatc <- function(gatc_sites.df = NULL, gatc.color = "red", gatc.size = 5,
    plot.space = 0.2, plot.height = 0.05) {
    structure(
        list(
            gatc_sites.df = gatc_sites.df, gatc.color = gatc.color, gatc.size = gatc.size,
            plot.space = plot.space, plot.height = plot.height
        ),
        class = "gatc"
    )
}


#' @export
ggplot_add.gatc <- function(object, plot, object_name) {
    if (!is.data.frame(object$gatc_sites.df)) {
        stop("data.frame of GATC sites is required")
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
    gatc_sites.df <- object$gatc_sites.df
    gatc.color <- object$gatc.color
    gatc.size <- object$gatc.size
    plot.space <- object$plot.space
    plot.height <- object$plot.height

    valid.bed <- ..getRegionsPlot(df = gatc_sites.df, columns = c("seqnames", "start", "end"), chr = plot.chr, start = plot.region.start, end = plot.region.end)

    gatc.plot <- ..plotPeak(valid.bed = valid.bed, plot.size = gatc.size, plot.color = gatc.color, peak.label = FALSE)

    gatc.plot <- gatc.plot +
        ggplot2::labs(y = "GATC") +
        ..peakGatcTheme(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end))

    # assemble plot
    patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
        gatc.plot,
        ncol = 1, heights = c(1, plot.height)
    )
}
