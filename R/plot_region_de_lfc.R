#' Plotting de results with lfc
#'
#' `geom_de.res.lfc` is a ggplot layer that visualises the de_results and logFC across a given region.
#' * regions are coloured by de result: 1, 0, -1, NA (grey for NA)
#' * cannot be plotted by itself, must be added to an existing plot - see examples.
#'
#'
#' @param de_results.df A data.frame of differential testing results as outputted from `edgeR_results()`
#' @param plot.space Specify gap to next plot, default is 0.1
#' @param plot.height Specify overall height of plot, default is 0.1
#'
#' @return ggplot_add object
#' @export
#'
#' @examples
#' counts.df <- random_counts()
#' dge <- edgeR_set_up(counts.df)
#' de_results <- random_edgeR_results()
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_de.res.lfc(de_results)
geom_de.res.lfc <- function(de_results.df,
                            plot.space = 0.1, plot.height = 0.1) {
  structure(list(
    de_results.df = de_results.df, plot.space = 0.1, plot.height = plot.height
  ),
  class = "de.res.lfc"
  )
}

#' @export
ggplot_add.de.res.lfc <- function(object, plot, object_name) {
  if(!is.data.frame(object$de_results.df)) {
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
  de_results.df <- object$de_results.df
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  df_regions <- dplyr::filter(de_results.df,
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

  df_fc <- df_regions %>%
    dplyr::summarise(abs_max = max(logFC),
                     abs_min = abs(min(logFC))) %>%
    dplyr::mutate(abs_fc = ifelse(abs_max >= abs_min, abs_max, abs_min),
                  abs_fc = round(abs_fc))

  colours <- c("Upreg" = "red", "Downreg" = "blue", "Not_included" = "grey", "No_sig" = "black")

  de_res.plot <- df_regions %>%
    ggplot2::ggplot(ggplot2::aes(x = start, y = de, colour = meth_status)) +
    ggplot2::geom_polygon(data = df_colour, ggplot2::aes(x = Position, y = y_axis_2, fill = meth_status)) +
    ggplot2::geom_segment(ggplot2::aes(xend=start, yend=0)) +
    ggplot2::geom_segment(ggplot2::aes(x=end, xend=end, y=de, yend=0)) +
    ggplot2::geom_segment(ggplot2::aes(x=start, xend=end, y=de, yend=de)) +
    ggplot2::geom_segment(ggplot2::aes(x=start, xend=end, y=0, yend=0)) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::coord_cartesian(xlim = c(plot.region.start, plot.region.end)) +
    ggplot2::scale_y_continuous(limits = c(-(df_fc$abs_fc) - 0.5, df_fc$abs_fc + 0.5),
                       expand = c(0, 0), position = "right") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(t = 0.1, b = 0.1)
    )

  # assemble plot
  patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                        de_res.plot,
                        ncol = 1, heights = c(1, plot.height)
  )
}
