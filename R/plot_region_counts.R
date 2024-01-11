#' Plot for counts for all samples across a given region
#'
#' `plot_counts_all_bams` plots a ggplot object visualising the raw counts from the bam files across a given region.
#' * this can be used as the base layer (set n_col = 1) for additional plot layers (geom_peak.new, geom_gatc, geom_de.res.lfc etc)
#'
#' @param counts.df A data.frame of counts as outputted from [process_bams()]
#' @param seqnames A character string of the chromosome of interest
#' @param start_region A number providing the start of region to plot
#' @param end_region A number providing the end of region to plot
#' @param n_col The number of columns to facet the graph by. Default is 1
#' @param layout Determines the layout of the plot. Default is "stacked" collapsing the Dam samples into one plot, and the Fusion samples into another. Samples can be plotted separately using "spread".
#' @param log2_scale Determines whether or not to display the counts on a log2 scale. Default is TRUE.
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
#'                      n_col = 1,
#'                      layout = "stacked")
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1,
#'                      layout = "spread")
#' # Can use this plot to layer other plots -----------------------------
#' dm_results <- random_edgeR_results()
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_dm.res.lfc(dm_results)
plot_counts_all_bams <- function(counts.df, seqnames, start_region = NULL, end_region = NULL, n_col = 1, layout = c("stacked", "spread"), log2_scale = TRUE) {
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
  df <- df %>% dplyr::filter(.data$seqnames == {{seqnames}}) %>%
    dplyr::filter(.data$start >= start_region, .data$end <= end_region)
  if(nrow(df) == 0) {
    stop("No data available for provided region, make the region larger")
  }
  df <- plot_counts_reshape(df)

  if(log2_scale == TRUE) {
    df$raw_counts <- log2(df$raw_counts)
  }
  if(length(layout) > 1) {
    layout <- "stacked"
  }
  if(layout == "stacked") {
    df <- df %>%
      dplyr::group_by(.data$dam) %>%
      dplyr::mutate(dam_num = as.integer(factor(.data$bam))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dam_2 = paste0(.data$dam, "_", .data$dam_num))

    plot <- counts_ggplot(df, start_region, end_region, "dam", seqnames, labs_fill = "Replicate", n_col=n_col, alpha = 0.5)
  } else if(layout == "spread") {

    plot <- counts_ggplot(df, start_region, end_region, "bam", seqnames, labs_fill = NULL, n_col=n_col) +
      ggplot2::scale_fill_discrete() +
      ggplot2::theme(legend.position = "none")
  }

  if(log2_scale == TRUE) {
    plot <- plot +
      ggplot2::labs(y = "Log2_counts")
  }
  plot
}

plot_counts_reshape <- function(counts) {
  df <- counts %>% dplyr::mutate(number = seq_len(dplyr::n())) %>%
    .[rep(seq_len(nrow(.)), times = 4),] %>%
    .[order(.$number),] %>%
    dplyr::group_by(.data$number) %>%
    dplyr::mutate(num = seq_len(dplyr::n())) %>%
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
                             colnames(.[, grepl(".bam", names(.))])) %>%
    dplyr::mutate(dam = ifelse(grepl("Dam", .data$bam), "Dam", "Fusion")) %>%
    .[order(.$dam, decreasing = TRUE),]
  df
}

counts_ggplot <- function(df, start_region, end_region, group, seqnames, labs_fill, n_col=n_col, ...) {
  if(!("dam_num" %in% colnames(df))) {
    df$dam_num <- 5
  }
  plot <- df %>%
    ggplot2::ggplot() +
    ggplot2::geom_polygon(ggplot2::aes(x = .data$Position, y = .data$raw_counts, fill = as.factor(.data$dam_num), ...), ...) +
    ggplot2::scale_fill_brewer() +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::coord_cartesian(xlim = c(start_region, end_region)) +
    ggplot2::facet_wrap(~ .data[[group]], ncol = n_col) +
    ggplot2::labs(title = paste0(seqnames, ":", start_region, "-", end_region), fill = labs_fill, ...) +
    ggplot2::theme_classic()
  plot
}



