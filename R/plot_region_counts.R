#' Plot for just counts for all samples across a given region
#'
#' @param df df of counts as outputted from [process_bams()]
#' @param seqnames chromosome of interest
#' @param start_region start
#' @param end_region end
#' @param n_col n_columns to facet the graph by - default is 3
#'
#' @return ggplot
#' @export
#'
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' plot_counts_all_bams(counts_df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 3)
#' # Can use this plot to layer other plots,
#'   just set n_col = 1 and start layering ----------------------------------
#' dge <- edgeR_set_up(counts.df)
#' de_results <- edgeR_results(dge, p.value = 0.05, lfc = 1)
#' de_results <- add_de(de_results, regions = regions_gatc_drosophila_dm6)
#' plot_counts_all_bams(counts_df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_de.res(de_results)
plot_counts_all_bams <- function(df, seqnames, start_region = NULL, end_region = NULL, n_col = 3) {
  if(!is.data.frame(df)) {
    stop("data.frame of counts is required")
  }
  if(!(seqnames %in% df$seqnames)) {
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
#' @return plot
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


#' @export
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
