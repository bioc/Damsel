#' Plot for just counts for all samples across a given region
#'
#' @param df df of counts as outputted from [process_bams()]
#' @param seqnames chromosome of interest
#' @param start_region start
#' @param end_region end
#' @param n_col n_columns to facet the graph by - default is 3
#'
#' @return
#' @export
#'
#' @examples
plot_counts_all_bams <- function(df, seqnames, start_region, end_region, n_col = 3) {
  df <- df
  colnames(df) <- chartr("-", "_", colnames(df))
  df <- df %>% dplyr::filter(seqnames == seqnames, start >= start_region, end <= end_region)
  df <- df %>% dplyr::mutate(number = 1:dplyr::n()) %>%
      .[rep(seq_len(nrow(.)), times = 4),] %>%
      .[order(.$number),] %>%
      dplyr::group_by(number) %>%
      dplyr::mutate(num = 1:dplyr::n()) %>%
      dplyr::mutate(Position = dplyr::case_when(num == 1 ~ start,
                                                num == 2 ~ start,
                                                num == 3 ~ end,
                                                TRUE ~ end))
  df <- df %>% dplyr::mutate_at(ggplot2::vars(matches(".bam")), ~ dplyr::case_when(num == 1 ~ 0,
                                                                                   num == 2 ~ .,
                                                                                   num == 3 ~ .,
                                                                                   TRUE ~ 0))
  df <- df %>% tidyr::gather(key = "bam",
                             value = "raw_counts",
                             colnames(.[, grepl(".bam", names(.))]))
  df %>%
    ggplot2::ggplot() +
    ggplot2::geom_polygon(ggplot2::aes(x = Position, y = raw_counts)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::facet_wrap(~ bam, ncol = n_col)
}


#' GGplot addable counts graph for all samples
#'
#' @param region.df df of counts
#' @param n_col n of columns to facet by
#' @param region.color colour
#' @param plot.space gap to next plot
#' @param plot.height height of plot
#'
#' @return
#' @export
#'
#' @examples
geom_regions.counts <- function(region.df = NULL, n_col = 1, region.color = "black",
                                plot.space = 0.1, plot.height = 1) {
  structure(list(
    region.df = region.df, n_col = n_col, region.color = region.color,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "regions.counts"
  )
}

#' Constructor for ggplot addable
#'
#' @param object regions.counts
#' @param plot plot being added to
#' @param object_name regions.counts
#'
#' @return
#' @export
#'
#' @examples
ggplot_add.regions.counts <- function(object, plot, object_name) {
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
  n_col <- object$n_col
  region.color <- object$region.color
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  df <- region.df
  colnames(df) <- chartr("-", "_", colnames(df))
  df <- df %>% dplyr::filter(seqnames == seqnames, start >= start_region, end <= end_region)
  df <- df %>% dplyr::mutate(number = 1:dplyr::n()) %>%
    .[rep(seq_len(nrow(.)), times = 4),] %>%
    .[order(.$number),] %>%
    dplyr::group_by(number) %>%
    dplyr::mutate(num = 1:dplyr::n()) %>%
    dplyr::mutate(Position = dplyr::case_when(num == 1 ~ start,
                                              num == 2 ~ start,
                                              num == 3 ~ end,
                                              TRUE ~ end))
  df <- df %>% dplyr::mutate_at(ggplot2::vars(matches(".bam")), ~ dplyr::case_when(num == 1 ~ 0,
                                                                                   num == 2 ~ .,
                                                                                   num == 3 ~ .,
                                                                                   TRUE ~ 0))

  counts.plot <- df %>%
    ggplot2::ggplot() +
    ggplot2::geom_polygon(ggplot2::aes(x = Position, y = raw_counts)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::facet_wrap(~ bam, ncol = n_col)

  patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                        counts.plot,
                        ncol = 1, heights = c(1, plot.height)
  )
}
