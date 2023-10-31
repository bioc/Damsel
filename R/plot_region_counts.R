#' Plot for counts for all samples across a given region
#'
#' `plot_counts_all_bams` plots a ggplot object visualising the raw counts from the bam files across a given region.
#' * this can be used as the base layer (set n_col = 1) for additional plot layers (geom_peak.new, geom_gatc, geom_de.res.lfc etc)
#'
#' @param counts.df A data.frame of counts as outputted from [process_bams()]
#' @param seqnames A character string of the chromosome of interest
#' @param start_region A number providing the start of region to plot
#' @param end_region A number providing the end of region to plot
#' @param n_col The number of columns to facet the graph by. Default is 3
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' set.seed(123)
#' counts.df <- random_counts()
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 3)
#' # Can use this plot to layer other plots -----------------------------
#' dm_results <- random_edgeR_results()
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_dm.res.lfc(dm_results)
plot_counts_all_bams <- function(counts.df, seqnames, start_region = NULL, end_region = NULL, n_col = 3) {
  if(!is.data.frame(counts.df)) {
    stop("data.frame of counts is required")
  }
  if(!(seqnames %in% counts.df$seqnames)) {
    stop("seqnames must be element of seqnames in provided data.frame")
  }
  if(is.null(start_region) | !is.numeric(start_region)) {
    stop("numeric element for start_region is required")
  }
  if(is.null(end_region) | !is.numeric(end_region)) {
    stop("numeric element for end_region is required")
  }
  if(end_region <= start_region) {
    stop("end_region must be greater than start_region")
  }
  df <- counts.df
  #colnames(df) <- chartr("-", "_", colnames(df))
  df <- df %>% dplyr::filter(.data$seqnames == {{seqnames}}) %>%
    dplyr::filter(.data$start >= start_region, .data$end <= end_region)
  if(nrow(df) == 0) {
    stop("No data available for provided region, make the region larger")
  }
  df <- df %>% dplyr::mutate(number = 1:dplyr::n()) %>%
      .[rep(seq_len(nrow(.)), times = 4),] %>%
      .[order(.$number),] %>%
      dplyr::group_by(.data$number) %>%
      dplyr::mutate(num = 1:dplyr::n()) %>%
      dplyr::mutate(Position = dplyr::case_when(.data$num == 1 ~ .data$start,
                                                .data$num == 2 ~ .data$start,
                                                .data$num == 3 ~ .data$end,
                                                TRUE ~ .data$end))
  df <- df %>% dplyr::mutate_at(ggplot2::vars(tidyr::matches(".bam")), ~ dplyr::case_when(.data$num == 1 ~ 0,
                                                                                   .data$num == 2 ~ .,
                                                                                   .data$num == 3 ~ .,
                                                                                   TRUE ~ 0))
  df <- df %>% tidyr::gather(key = "bam",
                             value = "raw_counts",
                             colnames(.[, grepl(".bam", names(.))]))
  df %>%
    ggplot2::ggplot() +
    ggplot2::geom_polygon(ggplot2::aes(x = .data$Position, y = .data$raw_counts)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::coord_cartesian(xlim = c(start_region, end_region)) +
    ggplot2::facet_wrap(~ bam, ncol = n_col) +
    ggplot2::labs(title = paste0(seqnames, ":", start_region, "-", end_region))
}



geom_regions.counts <- function(region.df = NULL, n_col = 1, region.color = "black",
                                plot.space = 0.1, plot.height = 1) {
  structure(list(
    region.df = region.df, n_col = n_col, region.color = region.color,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "regions.counts"
  )
}



ggplot_add.regions.counts <- function(object, plot, object_name) {
  # get plot data
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
  n_col <- object$n_col
  region.color <- object$region.color
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  df <- region.df
  colnames(df) <- chartr("-", "_", colnames(df))
  df <- df %>% dplyr::filter(.data$seqnames == plot.chr, .data$start >= plot.region.start, .data$end <= plot.region.end)
  df <- df %>% dplyr::mutate(number = 1:dplyr::n()) %>%
    .[rep(seq_len(nrow(.)), times = 4),] %>%
    .[order(.$number),] %>%
    dplyr::group_by(.data$number) %>%
    dplyr::mutate(num = 1:dplyr::n()) %>%
    dplyr::mutate(Position = dplyr::case_when(.data$num == 1 ~ .data$start,
                                              .data$num == 2 ~ .data$start,
                                              .data$num == 3 ~ .data$end,
                                              TRUE ~ .data$end))
  df <- df %>% dplyr::mutate_at(ggplot2::vars(tidyr::matches(".bam")), ~ dplyr::case_when(.data$num == 1 ~ 0,
                                                                                   .data$num == 2 ~ .,
                                                                                   .data$num == 3 ~ .,
                                                                                   TRUE ~ 0))
  counts.plot <- df %>%
    ggplot2::ggplot() +
    ggplot2::geom_polygon(ggplot2::aes(x = .data$Position, y = .data$raw_counts)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::coord_cartesian(xlim = c(plot.region.start, plot.region.end)) +
    ggplot2::facet_wrap(~ bam, ncol = n_col) +
    ggplot2::labs(title = paste0(seqnames, ":", plot.region.start, "-", plot.region.start))

  patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                        counts.plot,
                        ncol = 1, heights = c(1, plot.height)
  )
}
