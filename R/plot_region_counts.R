#' Plot for counts for all samples across a given region
#'
#' `plot_counts_all_bams` plots a ggplot2 object visualising the raw counts from the bam files across a given region.
#' * this can be used as the base layer (set n_col = 1) for additional plot layers (geom_peak.new, geom_gatc, geom_de.res.lfc etc)
#'
#' @param counts.df A data.frame of counts as outputted from [process_bams()].
#' @param seqnames A character string of the chromosome of interest.
#' @param start_region A number providing the start of region to plot.
#' @param end_region A number providing the end of region to plot.
#' @param layout Determines the layout of the plot. Default is "stacked" collapsing the Dam samples into one plot, and the Fusion samples into another. Samples can be plotted separately using "spread".
#' @param log2_scale Determines whether or not to display the counts on a log2 scale. Default is FALSE.
#' @param colours Specify colours for the replicates.
#' @param ... Arguments passed to ggplot2
#'
#' @return A `ggplot2` object.
#' @export
#' @references ggcoverage - Visualise and annotate omics coverage with ggplot2. https://github.com/showteeth/ggcoverage/tree/main
#' @seealso [geom_peak()] [geom_dm()] [geom_genes()] [geom_gatc()] [plot_wrap()]
#' @examples
#' set.seed(123)
#' counts.df <- random_counts()
#' plot_counts_all_bams(counts.df,
#'     seqnames = "chr2L",
#'     start_region = 1,
#'     end_region = 40000,
#'     layout = "stacked",
#'     log2_scale = FALSE
#' )
#' plot_counts_all_bams(counts.df,
#'     seqnames = "chr2L",
#'     start_region = 1,
#'     end_region = 40000,
#'     layout = "spread",
#'     log2_scale = FALSE
#' )
#' # Can use this plot to layer other plots -----------------------------
#' dm_results <- random_edgeR_results()
#' plot_counts_all_bams(counts.df,
#'     seqnames = "chr2L",
#'     start_region = 1,
#'     end_region = 40000,
#'     log2_scale = FALSE
#' ) +
#'     geom_dm.res.lfc(dm_results)
plot_counts_all_bams <- function(counts.df, seqnames, start_region=NULL, end_region=NULL, layout=c("stacked", "spread"), log2_scale=FALSE, colours=NULL, ...) {
    if (!is.data.frame(counts.df)) {
        stop("data.frame of counts is required")
    }
    if (!(seqnames %in% counts.df$seqnames)) {
        stop("seqnames must be element of seqnames in provided data.frame")
    }
    if (is.null(start_region) | !is.numeric(start_region)) {
        stop("numeric element for start_region is required")
    }
    if (is.null(end_region) | !is.numeric(end_region)) {
        stop("numeric element for end_region is required")
    }
    if (end_region <= start_region) {
        stop("end_region must be greater than start_region")
    }
    df <- counts.df
    df <- df %>%
        dplyr::filter(.data$seqnames == {{ seqnames }}) %>%
        dplyr::filter(.data$start >= start_region, .data$end <= end_region)
    if (nrow(df) == 0) {
        stop("No data available for provided region, make the region larger")
    }
    df <- ..plotCountsReshape(df)

    if (missing(colours)) {
        n_samples <- df %>%
            dplyr::ungroup() %>%
            dplyr::filter(.data$dam == "Dam") %>%
            dplyr::distinct(.data$bam) %>%
            nrow()
        if (n_samples == 2) {
            colours <- c("darkblue", "lightblue")
        } else if (n_samples == 3) {
            colours <- c("darkblue", "blue", "lightblue")
        }
    }
    if (log2_scale == TRUE) {
        df$raw_counts <- log2(df$raw_counts)
    }
    if (length(layout) > 1) {
        layout <- "stacked"
    }
    if (layout == "stacked") {
        df <- df %>%
            dplyr::group_by(.data$dam) %>%
            dplyr::mutate(dam_num = as.integer(factor(.data$bam))) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(dam_2 = paste0(.data$dam, "_", .data$dam_num))

        plot <- ..countsGgplot(df, start_region, end_region, colours, "dam",
            group_levels = c("Fusion", "Dam"),
            seqnames, labs_fill = "Replicate", alpha = 0.5
        )
    } else if (layout == "spread") {
        sample_order <- df$bam %>% unique()
        colours <- c("#480871", "#E30AD0")
        plot <- ..countsGgplot(df, start_region, end_region, colours, "bam",
            group_levels = sample_order,
            seqnames, labs_fill = "Replicate", alpha = 1
        ) +
            ggplot2::theme(legend.position = "none")
    }

    if (log2_scale == TRUE) {
        plot <- plot +
            ggplot2::labs(y = "Log2_counts")
    }
    plot
}

#' @export
#' @rdname plot_counts_all_bams
plotCounts <- plot_counts_all_bams

..plotCountsReshape <- function(counts) {
    df <- counts %>%
        dplyr::mutate(number = seq_len(dplyr::n())) %>%
        .[rep(seq_len(nrow(.)), times = 4), ] %>%
        .[order(.$number), ] %>%
        dplyr::group_by(.data$number) %>%
        dplyr::mutate(num = seq_len(dplyr::n())) %>%
        dplyr::mutate(Position = dplyr::case_when(
            .data$num == 1 ~ .data$start,
            .data$num == 2 ~ .data$start,
            .data$num == 3 ~ .data$end,
            TRUE ~ .data$end
        ))
    df <- df %>% dplyr::mutate_at(ggplot2::vars(tidyr::matches(".bam", ignore.case = TRUE)), ~ dplyr::case_when(
        .data$num == 1 ~ 0,
        .data$num == 2 ~ .,
        .data$num == 3 ~ .,
        TRUE ~ 0
    ))
    df <- df %>%
        tidyr::gather(
            key = "bam",
            value = "raw_counts",
            colnames(.[, grepl(".bam", names(.), ignore.case = TRUE)])
        ) %>%
        dplyr::mutate(dam = ifelse(grepl("Dam", .data$bam, ignore.case = TRUE), "Dam", "Fusion")) %>%
        .[order(.$dam, decreasing = TRUE), ]
    df
}

..countsGgplot <- function(df, start_region, end_region, colours, group, group_levels = NULL,
    seqnames, labs_fill, alpha, ...) {
    if (!("dam_num" %in% colnames(df))) {
        df$dam_num <- df$dam
    }
    plot <- df %>%
        ggplot2::ggplot() +
        ggplot2::geom_polygon(ggplot2::aes(x = .data$Position, y = .data$raw_counts, fill = as.factor(.data$dam_num)), alpha = alpha, ...) +
        ggplot2::scale_fill_manual(values = colours, ...) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(xlim = c(start_region, end_region)) +
        ggplot2::facet_wrap(~ factor(.data[[group]], levels = group_levels), ncol = 1) +
        ggplot2::labs(title = paste0(seqnames, ":", start_region, "-", end_region), fill = labs_fill, ...) +
        ggplot2::theme_classic()
    plot
}
