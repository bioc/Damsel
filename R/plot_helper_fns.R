## in plot_peak and plot_gatc

..getRegionsPlot <- function(df, columns, chr, start, end = NULL) {
    # subset used chromosome
    df <- df[, columns]
    df$start <- as.numeric(df$start) + 1
    df <- df[df$seqnames == chr, ] %>% dplyr::arrange(start)
    rownames(df) <- NULL

    df.select <- df[df$end >= start & df$start <= end, ]
    if (nrow(df.select) == 0) {
        return(df.select)
    }
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

..plotPeak <- function(valid.bed, plot.size, plot.color, peak.label) {
    if (nrow(valid.bed) == 0) {
        message("No data available for this region")
        peak.plot <- ggplot2::ggplot() +
            ggplot2::geom_blank()
    } else {
        peak.plot <- valid.bed %>%
            ggplot2::ggplot(ggplot2::aes(x = .data$start, y = 1)) +
            ggplot2::geom_segment(
                mapping = ggplot2::aes(
                    x = .data$start, y = 1,
                    xend = .data$end, yend = 1
                ),
                linewidth = plot.size, color = plot.color
            ) +
            if (peak.label == TRUE) {
                ggplot2::geom_label(ggplot2::aes(
                    x = (.data$start + .data$end) / 2,
                    label = .data$peak_id
                ), colour = "black", size = 3, vjust = "bottom", nudge_y = 0.02)
            }
    }
    peak.plot
}

..peakGatcTheme <- function(margin.len, x.range) {
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
