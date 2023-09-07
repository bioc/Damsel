#' Plot for a GATC track
#'
#' @param gatc.df df of positions of GATC sites - can be made from GATC region file
#' @param gatc.color default is red
#' @param gatc.size size of line
#' @param plot.space gap to next plot
#' @param plot.height size of plot
#'
#' @return plot
#' @export
#'
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' gatc_sites <- dplyr::mutate(regions_gatc_drosophila_dm6,
#'                             seqnames = paste0("chr", seqnames),
#'                             start = start - 3, end = start + 4)
#'
#' plot_counts_all_bams(counts_df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_gatc(gatc_sites)
#' # The plots can be layered -------------------------------------------------
geom_gatc <- function(gatc.df = NULL, gatc.color = "red", gatc.size = 5,
                      plot.space = 0.1, plot.height = 0.1) {
  structure(list(
    gatc.df = gatc.df, gatc.color = gatc.color, gatc.size = gatc.size,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "gatc"
  )
}


#' @export
ggplot_add.gatc <- function(object, plot, object_name) {
  if(!is.data.frame(object$gatc.df)) {
    stop("data.frame of GATC sites is required")
  }
  # get plot data
  # get plot data, plot data should contain bins
  if (("patchwork" %in% class(plot)) && length(plot[[1]]$layers) != 2) {
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
  gatc.df <- object$gatc.df
  gatc.color <- object$gatc.color
  gatc.size <- object$gatc.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  bed.info <- gatc.df
  bed.info <- bed.info[,c("seqnames", "start", "end")]
  # convert to 1-based
  bed.info$start <- as.numeric(bed.info$start) + 1

  # get valid bed
  valid.bed <- GetRegion_hack(chr = plot.chr, df = bed.info, start = plot.region.start, end = plot.region.end)

  gatc.plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = valid.bed,
      mapping = ggplot2::aes_string(
        x = "start",
        y = "1",
        xend = "end",
        yend = "1"
      ),
      size = gatc.size,
      color = gatc.color
    ) +
    ggplot2::labs(y = "GATC")

  # add theme
  gatc.plot <- gatc.plot + theme_peak_hack(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end))

  # assemble plot
  patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                        gatc.plot,
                        ncol = 1, heights = c(1, plot.height)
  )
}


