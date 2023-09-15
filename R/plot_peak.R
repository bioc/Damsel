#' Plotting peaks
#'
#' `geom_peak.new` is a ggplot layer that visualises the positions of peaks across a given region.
#' * cannot be plotted by itself, must be added to an existing ggplot object - see examples.
#'
#'
#' @param peak.df data frame of peak locations
#' @param peak.color can set the colour if you want
#' @param peak.size default is 5
#' @param plot.space gap to next plot - default 0.1
#' @param plot.height height of plot - leave it at default
#'
#' @return ggplot_add object
#' @export
#'
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' dge <- edgeR_set_up(counts.df)
#' de_results <- edgeR_results(dge, p.value = 0.05, lfc = 1)
#' peaks <- aggregate_peaks(de_results, regions = regions_gatc_drosophila_dm6)
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_peak.new(peaks)
#' # The plots can be layered -------------------------------------------------
geom_peak.new <- function(peak.df = NULL, peak.color = "black", peak.size = 5,
                          plot.space = 0.1, plot.height = 0.05) {
  structure(list(
    peak.df = peak.df, peak.color = peak.color, peak.size = peak.size,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "peak.new"
  )
}


#' @export
ggplot_add.peak.new <- function(object, plot, object_name) {
  if(!is.data.frame(object$peak.df)) {
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
  bed.file <- object$bed.file
  peak.df <- object$peak.df
  peak.color <- object$peak.color
  peak.size <- object$peak.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  bed.info <- peak.df
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
    patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
                          peak.plot,
                          ncol = 1, heights = c(1, plot.height))
  }

  peak.plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = valid.bed,
      mapping = ggplot2::aes(
        x = start,
        y = 1,
        xend = end,
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

