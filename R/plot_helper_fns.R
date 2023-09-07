##in plot_peak and plot_gatc

GetRegion_hack <- function(df, chr, start, end = NULL) {
  # subset used chromosome
  df <- df[df$seqnames == chr, ] %>% dplyr::arrange(start)
  rownames(df) <- NULL

  df.select <- df[df$end >= start & df$start <= end, ]
  init.start <- df.select[1, "start"]
  if (init.start < start) {
    df.select[1, "start"] <- start
  }
  if (!is.null(end)) {
    final.end <- df.select[nrow(df.select), "end"]
    if (final.end > end) {
      df.select[nrow(df.select), "end"] <- end
    }
  }
  return(df.select)
}


theme_peak_hack <- function(margin.len, x.range) {
  list(
    ggplot2::theme_classic(),
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(t = margin.len, b = margin.len)
    ),
    ggplot2::scale_y_continuous(
      limits = c(1 - 0.1, 1 + 0.1),
      expand = c(0, 0), position = "right"
    ),
    ggplot2::scale_x_continuous(expand = c(0, 0)),
    ggplot2::coord_cartesian(xlim = x.range)
    )
}
