#' Plotting peaks
#'
#' `geom_peak.new` adds the peak location to an existing ggplot2 object
#'
#' @param bed.file bed input of peaks
#' @param peak.df data frame of peak locations
#' @param peak.color can set the colour if you want
#' @param peak.size default is 5
#' @param plot.space gap to next plot - default 0.1
#' @param plot.height height of plot - leave it at default
#'
#' @return
#' @export
#'
#' @examples
geom_peak.new <- function(bed.file = NULL, peak.df = NULL, peak.color = "black", peak.size = 5,
                          plot.space = 0.1, plot.height = 0.1) {
  structure(list(
    bed.file = bed.file, peak.df = peak.df, peak.color = peak.color, peak.size = peak.size,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "peak.new"
  )
}

#' how to make the above plot
#'
#' code for the `geom_peak.new` plot
#'
#' @param object this is what goes in - peak.df or bed.file
#' @param plot this is the plot you're adding to
#' @param object_name "peak.new
#'
#' @return
#' @export
#'
#' @examples
ggplot_add.peak.new <- function(object, plot, object_name) {
  # get plot data
  # get plot data, plot data should contain bins
  if (("patchwork" %in% class(plot)) && length(plot[[1]]$layers) == 1) {
    plot.data <- plot[[1]]$layers[[1]]$data
    if(!("data.frame" %in% class(plot.data))) {
      plot.data <- plot[[1]]$data
    }
  } else if ("patchwork" %in% class(plot) && length(plot[[1]]$layers) == 2) {
    plot.data <- plot[[1]]$layers[[2]]$data
    colnames(plot.data) <- c("start", "end", "y1", "y2", "seqnames")
  } else if (!("patchwork" %in% class(plot)) && length(plot$layers) == 1) {
    plot.data <- plot$layers[[1]]$data
    if(!("data.frame" %in% class(plot.data))) {
      plot.data <- plot$data
    }
  } else if (!("patchwork" %in% class(plot)) && length(plot$layers) == 2) {
    plot.data <- plot$layers[[2]]$data
    colnames(plot.data) <- c("start", "end", "y1", "y2", "seqnames")
  }
  # prepare plot range
  # the plot region are not normal, so start is minimum value
  plot.chr <- as.character(plot.data[1, "seqnames"])
  plot.region.start <- min(plot.data[, "start"])
  plot.region.end <- max(plot.data[, "end"])

  # get parameters
  bed.file <- object$bed.file
  peak.df <- object$peak.df
  peak.color <- object$peak.color
  peak.size <- object$peak.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # prepare peak dataframe
  if (!is.null(bed.file)) {
    bed.info <- utils::read.table(file = bed.file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    bed.info <- bed.info[c(1, 2, 3)]
    colnames(bed.info) <- c("seqnames", "start", "end")
  } else if (!is.null(peak.df)) {
    bed.info <- peak.df
    bed.info <- bed.info[,c("seqnames", "start", "end")]
  }
  # convert to 1-based
  bed.info$start <- as.numeric(bed.info$start) + 1

  # get valid bed
  valid.bed <- GetRegion_hack(chr = plot.chr, df = bed.info, start = plot.region.start, end = plot.region.end)

  peak.plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = valid.bed,
      mapping = ggplot2::aes_string(
        x = "start",
        y = "1",
        xend = "end",
        yend = "1"
      ),
      size = peak.size,
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

