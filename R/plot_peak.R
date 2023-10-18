#' Plotting peaks
#'
#' `geom_peak.new` is a ggplot layer that visualises the positions of peaks across a given region.
#' * cannot be plotted by itself, must be added to an existing ggplot object - see examples.
#'
#'
#' @param peak.df A data.frame of peaks as outputted from `aggregate_peaks()`
#' @param peak.color Specify colour of peak, default is black
#' @param peak.size Specify size of rectangle, default is 5
#' @param plot.space Specify gap to next plot, default is 0.1
#' @param plot.height Specify overall height of plot, default is 0.05
#'
#' @return ggplot_add object
#' @export
#'
#' @examples
#' set.seed(123)
#' counts.df <- random_counts()
#' dm_results <- random_edgeR_results()
#' peaks <- aggregate_peaks(dm_results)
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_peak.new(peaks)
#' # The plots can be layered -------------------------------------------------
geom_peak.new <- function(peaks.df = NULL, peak.color = "black", peak.size = 5,
                          plot.space = 0.1, plot.height = 0.05) {
  structure(list(
    peaks.df = peaks.df, peak.color = peak.color, peak.size = peak.size,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "peak.new"
  )
}


#' @export
ggplot_add.peak.new <- function(object, plot, object_name) {
  if(!is.data.frame(object$peaks.df)) {
    stop("data.frame of peaks is required")
  }
  plot2 <- plot
  while("patchwork" %in% class(plot2)) {
    plot2 <- plot2[[1]]
  }
  plot.data <- plot2$labels$title
  plot.data <- stringr::str_split_1(plot.data, ":")
  # prepare plot range
  # the plot region are not normal, so start is minimum value
  plot.chr <- plot.data[1]
  plot.data <- stringr::str_split_1(plot.data[2], "-")
  plot.region.start <- plot.data[1] %>% as.numeric()
  plot.region.end <- plot.data[2] %>% as.numeric()

  # get parameters
  peaks.df <- object$peaks.df
  peak.color <- object$peak.color
  peak.size <- object$peak.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  bed.info <- peaks.df
  bed.info <- bed.info[,c("seqnames", "start", "end")]

  # convert to 1-based
  bed.info$start <- as.numeric(bed.info$start) + 1

  # get valid bed
  valid.bed <- GetRegion_hack(chr = plot.chr, df = bed.info, start = plot.region.start, end = plot.region.end)
  if(nrow(valid.bed) == 0) {
    message("No peak data available for this region")
    peak.plot <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::labs(y = "Peak") +
      theme_peak_hack(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end))
    patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                          peak.plot,
                          ncol = 1, heights = c(1, plot.height))
  }

  peak.plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = valid.bed,
      mapping = ggplot2::aes(
        x = .data$start,
        y = 1,
        xend = .data$end,
        yend = 1
      ),
      linewidth = peak.size,
      color = peak.color
    ) +
    ggplot2::labs(y = "Peak")

  # add theme
  peak.plot <- peak.plot + theme_peak_hack(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end))
  # assemble plot
  patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                        peak.plot,
                        ncol = 1, heights = c(1, plot.height)
  )
}

