geom_regions.lfc <- function(region.df = NULL, region.color = "black",
                             plot.space = 0.1, plot.height = 1) {
  structure(list(
    region.df = region.df, region.color = region.color,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "regions.lfc"
  )
}

ggplot_add.regions.lfc <- function(object, plot, object_name) {
  # get plot data
  # get plot data, plot data should contain bins
  if (("patchwork" %in% class(plot)) && length(plot[[1]]$layers) == 1) {
    plot.data <- plot[[1]]$layers[[1]]$data
  } else if ("patchwork" %in% class(plot) && length(plot[[1]]$layers) == 2) {
    plot.data <- plot[[1]]$layers[[2]]$data
    colnames(plot.data) <- c("start", "end", "y1", "y2", "seqnames")
  } else if (!("patchwork" %in% class(plot)) && length(plot$layers) == 1) {
    plot.data <- plot$layers[[1]]$data
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
  region.df <- object$region.df
  region.color <- object$region.color
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  #df_regions <- dd_de(final_edgeR_results_ad, regions = regions_for_testing)
  df_regions <- dplyr::filter(region.df, seqnames == plot.chr, start >= plot.region.start, end <= plot.region.end)
  df_colour <- df_regions %>%
      dplyr::mutate(number = 1:dplyr::n()) %>%
      .[rep(seq_len(nrow(.)), times = 4),] %>%
      .[order(.$number),] %>%
      dplyr::group_by(number) %>%
      dplyr::mutate(num = 1:dplyr::n()) %>%
      dplyr::mutate(Position = dplyr::case_when(num == 1 ~ start,
                                  num == 2 ~ start,
                                  num == 3 ~ end,
                                  TRUE ~ end),
             y_axis_2 = dplyr::case_when(num == 1 ~ 0,
                                  num == 2 ~ logFC,
                                  num == 3 ~ logFC,
                                  TRUE ~ 0))

  regions.plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = df_colour, ggplot2::aes(x = Position, y = y_axis_2)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::labs(y = "logFC")
  patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                        regions.plot,
                        ncol = 1, heights = c(1, plot.height)
  )
}
