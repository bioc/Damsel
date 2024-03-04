#' Identify peaks from differentially methylated regions
#'
#' `identifyPeaks` aggregates neighbouring differentially methylated regions, identifying 'peaks' where the provided transcription factor is believed to have bound to the DNA. These locations can then be used to identify the potential target genes.
#'
#' Small unmethylated regions are able to be 'skipped' over and included into peaks through the gap_size paramter, whose default is 150bp. This was selected due to the common approach of 75bp sequencing of DamID from the edges of the fragments.
#' The FDR and logFC for each peak is calculated via the theory of [csaw::getBestTest()] where the 'best' (smallest) p-value in the regions that make up the peak is selected as representative of the peak. The logFC is therefore the corresponding logFC from the FDR.
#'
#' @param dm_results The results from differential testing.
#' @param gap_size The maximum gap in base pairs between differentially methylated regions to be 'skipped'. Default is 150
#'
#' @return A data.frame of peaks ranked by p-value.
#' @export
#' @references Lun ATL, Smyth GK (2016). “csaw: a Bioconductor package for differential binding analysis of ChIP-seq data using sliding windows.” Nucleic Acids Res., 44(5), e45.
#' Lun ATL, Smyth GK (2014). “De novo detection of differentially bound regions for ChIP-seq data using peaks and windows: controlling error rates correctly.” Nucleic Acids Res., 42(11), e95.
#' @examples
#' set.seed(123)
#' counts.df <- random_counts()
#' dm_results <- random_edgeR_results()
#' peaks <- identifyPeaks(dm_results)
#' peaks
identifyPeaks <- function(dm_results, gap_size=150) {
    if (!is.data.frame(dm_results)) {
        stop("Must have data frame of differential_testing results from `edgeR_results`")
    }
    df <- data.frame(dm_results)
    peaks <- df %>%
        dplyr::mutate(
            number = seq_len(nrow(.)),
            a = unsplit(
                lapply(
                    split(.[, "meth_status"], .data$seqnames),
                    function(x) {
                        sequence(rle(x)$lengths)
                    }
                ),
                .data$seqnames
            ),
            b = ifelse(dplyr::lead(.data$a) == .data$a, 0, .data$a),
            multiple = ifelse(.data$b == 0, FALSE, TRUE)
        ) %>%
        dplyr::filter(.data$dm == 1) %>%
        dplyr::mutate(
            gap_regions = ifelse(dplyr::lead(.data$number) != .data$number + 1,
                dplyr::lead(.data$number) - .data$number - 1, NA
            ),
            gap_width = ifelse(is.na(.data$gap_regions), NA,
                dplyr::lead(.data$start) - .data$end
            ),
            gap_regions = replace(.data$gap_regions, dplyr::n(), 0),
            gap_width = replace(.data$gap_width, dplyr::n(), 0)
        ) %>%
        dplyr::mutate(
            new_one = ifelse(dplyr::lead(.data$a) < .data$a, "break", NA),
            new_one = ifelse(dplyr::lead(.data$a) == .data$a, "single", .data$new_one),
            new_one = replace(.data$new_one, dplyr::n(), "break")
        ) %>%
        dplyr::group_by(.data$new_one) %>%
        dplyr::mutate(seq = ifelse(.data$new_one == "break", seq_len(nrow(.)), 0)) %>%
        dplyr::ungroup() %>%
        tidyr::fill(seq, .direction = "up")

    peaks <- ..simplifyPeaks(peaks)
    # gaps <- gaps_new(df = peaks, dm_results = df, gap_size = {{ gap_size }})

    gaps <- ..checkForGaps(df = peaks, gap_size = {{ gap_size }})
    if (nrow(gaps) != 0) {
        #  gaps <- ..peakGaps(df = peaks, dm_results = df, gap_size = {{ gap_size }})
        gaps_ <- ..findGaps(gaps, gap_size = {{ gap_size }})
        gaps_ <- ..simplifyGaps(gaps_, gaps)
        gaps <- ..overlapGaps(gaps_, peaks, dm_results)
    }
    combo <- ..peaksCombine(peaks, gaps)
    combo <- data.frame(combo)
    combo
}


..simplifyPeaks <- function(peaks) {
    df <- peaks
    df <- df %>%
        dplyr::group_by(.data$seq) %>%
        dplyr::mutate(
            start = ifelse(.data$seq != 0, min(.data$start), .data$start),
            end = ifelse(.data$seq != 0, max(.data$end), .data$end),
            width = .data$end - .data$start + 1,
            FDR = ifelse(.data$seq != 0, min(.data$adjust.p), .data$adjust.p),
            n_regions_dm = ifelse(.data$seq != 0, dplyr::n(), 1),
            logFC_match = ifelse(.data$adjust.p == .data$FDR, .data$logFC, NA),
            region_pos = ifelse(.data$adjust.p == .data$FDR, .data$Position, NA)
        ) %>%
        tidyr::fill(c(logFC_match, region_pos, gap_regions, gap_width), .direction = "downup")

    peaks_num <- df %>%
        dplyr::ungroup() %>%
        dplyr::summarise(max_seq = max(.data$seq, na.rm = TRUE)) %>%
        .$max_seq

    df <- df %>%
        dplyr::mutate(seq = ifelse(.data$seq == 0,
            seq(from = peaks_num + 1, to = dplyr::n()), .data$seq
        )) %>%
        dplyr::filter(dplyr::row_number() == 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            peak_id = seq_len(nrow(.)),
        ) %>%
        .[, c(
            "peak_id", "Position", "seqnames", "start", "end", "width",
            "strand", "logFC_match", "FDR",
            "n_regions_dm", "region_pos", "gap_regions", "gap_width"
        )]
    df
}


..checkForGaps <- function(df, gap_size) {
    gaps <- df %>%
        dplyr::mutate(
            gap_width = dplyr::coalesce(.data$gap_width, 100000),
            gap_width = ifelse(.data$gap_width < 0, 10000, .data$gap_width)
        ) %>%
        dplyr::group_by(.data$seqnames) %>%
        dplyr::mutate(total_peaks = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(gap_width = ifelse(.data$total_peaks == 1, -5, .data$gap_width)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$gap_width > 0) %>%
        .[, !colnames(.) == "total_peaks"]
    if (nrow(gaps) == 0) {
        gaps$multiple_peaks <- numeric(0)
        gaps$n_regions_not_dm <- numeric(0)
        gaps$id_list <- character(0)
        gaps <- gaps[, c(
            "peak_id", "Position", "seqnames", "start",
            "end", "width", "strand", "logFC_match", "FDR",
            "multiple_peaks", "n_regions_dm",
            "n_regions_not_dm", "region_pos", "id_list"
        )]
        return(gaps)
    }
    gaps_check <- gaps %>%
        dplyr::mutate(gap2 = ifelse(.data$gap_width <= {{ gap_size }} & .data$width <= 10000, TRUE, FALSE)) %>%
        dplyr::filter(.data$gap2 != FALSE) %>%
        .[, !colnames(.) == "gap2"]

    if (nrow(gaps_check) == 0) {
        gaps_check$multiple_peaks <- numeric(0)
        gaps_check$n_regions_not_dm <- numeric(0)
        gaps_check$id_list <- character(0)
        gaps_check <- gaps_check[, c(
            "peak_id", "Position", "seqnames", "start",
            "end", "width", "strand", "logFC_match", "FDR",
            "multiple_peaks", "n_regions_dm",
            "n_regions_not_dm", "region_pos", "id_list"
        )]

        return(gaps_check)
    }
    gaps <- data.frame(gaps)
    gaps
}

..findGaps <- function(gaps, gap_size) {
    gaps_work <- rbind(c(1, 1, 1, 1)) %>%
        data.frame() %>%
        stats::setNames(c("seqnames", "id_list", "start", "end")) %>%
        dplyr::mutate(seqnames = as.character(.data$seqnames))
    end <- c()
    i <- 1
    number <- c(seq_len(nrow(gaps)))
    j <- i
    for (i in number) {
        if (gaps[i, ]$gap_width <= gap_size && gaps[i, ]$width < 10000) {
            start <- gaps[i, ]$start
            id_list <- gaps[i, ]$peak_id
            n_regions_dm <- gaps[i, ]$n_regions_dm
            n_regions_not_dm <- gaps[i, ]$gap_regions
            for (j in (i + 1):nrow(gaps)) {
                end <- c(end, gaps[j, ]$end)
                id_list <- paste0(id_list, ",", gaps[j, ]$peak_id)
                n_regions_dm <- n_regions_dm + gaps[j, ]$n_regions_dm
                n_regions_not_dm <- n_regions_not_dm + gaps[j, ]$gap_regions
                if (gaps[j, ]$gap_width > 150) {
                    break
                }
            }
            gaps_work[nrow(gaps_work) + 1, ] <- c(gaps[utils::tail(i, 1), ]$seqnames, id_list, start, utils::tail(end, 1))
            if (utils::tail(number, 1) == 100) {
                break
            }
        }
    }
    gaps_work <- data.frame(gaps_work)
    gaps_work
}

..simplifyGaps <- function(gaps_work, gaps) {
    gaps_work <- gaps_work[2:nrow(gaps_work), ] %>%
        dplyr::group_by(.data$end) %>%
        dplyr::mutate(
            number = seq_len(dplyr::n()),
            start = as.numeric(.data$start),
            end = as.numeric(.data$end)
        ) %>%
        dplyr::filter(.data$number == 1) %>%
        data.frame() %>%
        .[, colnames(.) != "number"] %>%
        dplyr::mutate(n_peaks = stringr::str_count(.data$id_list, ",") + 1)
    gaps_work <- gaps_work %>%
        dplyr::mutate(
            first = substr(.data$id_list, 1, regexpr(",", .data$id_list) - 1),
            first = as.double(.data$first)
        ) %>%
        data.frame()
    gaps_work$seqnames <- gaps[match(gaps_work$first, gaps$peak_id), "seqnames"]
    gaps_work$width <- gaps_work$end - gaps_work$start + 1
    gaps_work <- data.frame(gaps_work)
    gaps_work
}

..overlapGaps <- function(gaps_work, df, dm_results) {
    overlap_gaps_peaks <- plyranges::find_overlaps(plyranges::as_granges(gaps_work), plyranges::as_granges(df)) %>%
        data.frame() %>%
        dplyr::mutate(
            peak_id = .data$first, multiple_peaks = .data$n_peaks,
            Position = paste0(.data$seqnames, "-", .data$start),
            peak_id = paste0("PM_", .data$peak_id)
        ) %>%
        .[, c(
            "peak_id", "Position", "seqnames", "start", "end", "width", "strand",
            "logFC_match", "FDR", "multiple_peaks", "region_pos", "id_list"
        )] %>%
        dplyr::group_by(.data$Position) %>%
        dplyr::filter(.data$FDR == min(.data$FDR)) %>%
        dplyr::ungroup()

    overlap_gaps <- plyranges::find_overlaps(plyranges::as_granges(overlap_gaps_peaks), plyranges::as_granges(dm_results)) %>%
        data.frame() %>%
        dplyr::mutate(is_dm = ifelse(.data$dm == 1, "Dm", NA)) %>%
        dplyr::group_by(.data$peak_id, .data$is_dm) %>%
        dplyr::mutate(
            n_regions_dm = dplyr::n(), n_regions_not_dm = dplyr::n(),
            n_regions_dm = ifelse(is.na(.data$is_dm), NA, .data$n_regions_dm),
            n_regions_not_dm = ifelse(is.na(.data$is_dm), .data$n_regions_not_dm, NA)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$peak_id) %>%
        tidyr::fill(n_regions_dm, .direction = "downup") %>%
        tidyr::fill(n_regions_not_dm, .direction = "updown") %>%
        dplyr::filter(dplyr::row_number() == 1) %>%
        dplyr::ungroup()
    names(overlap_gaps) <- names(overlap_gaps) %>% gsub(".x", "", .)
    overlap_gaps <- overlap_gaps[, c(
        "peak_id", "Position", "seqnames", "start",
        "end", "width", "strand", "logFC_match", "FDR",
        "multiple_peaks", "n_regions_dm",
        "n_regions_not_dm", "region_pos", "id_list"
    )]

    overlap_gaps <- data.frame(overlap_gaps)
    overlap_gaps
}

..peaksCombine <- function(peaks, gaps) {
    peaks.df <- peaks
    gaps.df <- gaps

    peaks.df <- peaks.df %>%
        dplyr::filter(!.data$peak_id %in% as.numeric(unlist(strsplit(as.character(gaps.df$id_list), ","))))

    peaks.df <- peaks.df %>%
        dplyr::mutate(
            multiple_peaks = 1,
            n_regions_not_dm = 0,
            peak_id = paste0("PS_", .data$peak_id)
        ) %>%
        .[, c(
            "peak_id", "seqnames", "start", "end", "width",
            "strand", "logFC_match", "FDR",
            "multiple_peaks", "region_pos", "n_regions_dm", "n_regions_not_dm"
        )]
    gaps.df <- gaps.df[, c(
        "peak_id", "seqnames", "start", "end", "width",
        "strand", "logFC_match", "FDR", "multiple_peaks", "region_pos",
        "n_regions_dm", "n_regions_not_dm"
    )]

    df <- rbind(peaks.df, gaps.df)
    df <- df %>%
        dplyr::filter(.data$n_regions_dm > 2) %>%
        .[order(.$logFC_match, decreasing = TRUE), ] %>%
        .[order(.$FDR), ] %>%
        dplyr::mutate(rank_p = seq_len(nrow(.)))

    df <- df[, c(
        "peak_id", "seqnames", "start", "end", "width", "strand",
        "rank_p", "logFC_match", "FDR", "multiple_peaks", "region_pos",
        "n_regions_dm", "n_regions_not_dm"
    )]
    df
}
