#' Plot lfc results across the regions
#'
#' `geom_regions.lfc` is a ggplot layer that visualises the logFC results across a given region.
#' * cannot be plotted by itself, must be added to an existing plot - see examples.
#' * recommend using geom_de.res.lfc instead
#'
#'
#' @param region.df results df plus add_de: `add_de(edgeR_results())`
#' @param region.color leave as default
#' @param plot.space gap to next plot
#' @param plot.height height of plot
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
#' de_results <- add_de(de_results, regions = regions_gatc_drosophila_dm6)
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_regions.lfc(de_results)
#' # The plots can be layered -------------------------------------------------
geom_regions.lfc <- function(region.df, region.color = "black",
                             plot.space = 0.1, plot.height = 0.2) {
  structure(list(
    region.df = region.df, region.color = region.color,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "regions.lfc"
  )
}


#' @export
ggplot_add.regions.lfc <- function(object, plot, object_name) {
  if(!is.data.frame(object$region.df)) {
    stop("data.frame of de results is required")
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
  region.df <- object$region.df
  region.color <- object$region.color
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  #df_regions <- dd_de(final_edgeR_results_ad, regions = regions_for_testing)
  df_regions <- dplyr::filter(region.df,
                              seqnames == plot.chr,
                              start >= plot.region.start,
                              end <= plot.region.end)
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
