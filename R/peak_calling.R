new_fn_peaks <- function(dm_results, gap_size = 150) {
    df <- dm_results
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
            )
        ) %>%
        dplyr::mutate(
            new_one = ifelse(dplyr::lead(.data$a) < .data$a, "break", NA),
            new_one = ifelse(dplyr::lead(.data$a) == .data$a, "single", .data$new_one)
        ) %>%
        dplyr::group_by(.data$new_one) %>%
        dplyr::mutate(seq = ifelse(.data$new_one == "break", seq_len(nrow(.)), 0)) %>%
        dplyr::ungroup() %>%
        tidyr::fill(seq, .direction = "up")

    peaks <- simplify_peaks(peaks)

    gaps <- new_fn_gaps(peaks, {{ gap_size }})

    combo <- new_fn_p_combine(peaks, gaps)
    combo <- data.frame(combo)
    combo
}

simplify_peaks <- function(peaks) {
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
            multiple_peaks = NA
        ) %>%
        .[, c(
            "peak_id", "Position", "seqnames", "start", "end", "width",
            "strand", "logFC_match", "FDR", "multiple_peaks",
            "n_regions_dm", "region_pos", "gap_regions", "gap_width"
        )]
    df
}

new_fn_gaps <- function(peaks, gap_size = 150) {
    df <- peaks
    df <- df %>%
        dplyr::mutate(
            gap_width = ifelse(.data$gap_width < 0, 1e6, .data$gap_width),
            less_gap = ifelse(.data$gap_width <= gap_size, "Yes", "No"),
            a = sequence(rle(.data$less_gap)$lengths),
            b = ifelse(dplyr::lead(.data$a) == .data$a, 0, .data$a),
            multiple = ifelse(.data$b == 0, FALSE, TRUE)
        ) %>%
        dplyr::mutate(
            end = ifelse(.data$multiple == FALSE & .data$less_gap == "Yes",
                dplyr::lead(.data$end), NA
            ),
            n_regions_dm = ifelse(.data$multiple == FALSE & .data$less_gap == "Yes",
            .data$n_regions_dm + dplyr::lead(.data$n_regions_dm), NA
            ),
            n_regions_not_dm = ifelse(.data$multiple == FALSE &
                .data$less_gap == "Yes", .data$gap_regions, NA
            ),
            id_list = ifelse(.data$multiple == FALSE & .data$less_gap == "Yes",
                paste0(.data$peak_id, ",", dplyr::lead(.data$peak_id)), NA
            ),
            id_list = as.character(id_list)
        ) %>%
        dplyr::mutate(
            min_p_combo = ifelse(.data$multiple == FALSE & .data$less_gap == "Yes",
                dplyr::lead(.data$FDR), NA
            ),
            min_p_combo = pmin(.data$FDR, .data$min_p_combo),
            FDR = ifelse(.data$FDR == .data$min_p_combo, .data$FDR,
                .data$min_p_combo
            ),
            logFC_match = ifelse(.data$FDR == .data$min_p_combo, .data$logFC_match,
                dplyr::lead(.data$logFC_match)
            ),
            region_pos = ifelse(.data$FDR == .data$min_p_combo,
                .data$region_pos, dplyr::lead(.data$region_pos)
            )
        ) %>%
        dplyr::filter(!is.na(.data$end)) %>%
        dplyr::mutate(
            width = .data$end - .data$start + 1,
            peak_id = paste0("PM_", .data$peak_id),
            multiple_peaks = lengths(strsplit(id_list, ","))
        ) %>%
        .[, c(
            "peak_id", "Position", "seqnames", "start", "end", "width",
            "strand", "logFC_match", "FDR", "multiple_peaks",
            "n_regions_dm", "n_regions_not_dm", "region_pos", "id_list"
        )]
    df
}

new_fn_p_combine <- function(peaks, gaps) {
    peaks.df <- peaks
    gaps.df <- gaps

    peaks.df <- peaks.df %>%
        dplyr::filter(!.data$peak_id %in% unlist(strsplit(gaps.df$id_list, ",")))

    peaks.df <- peaks.df[, c(
        "peak_id", "seqnames", "start", "end", "width",
        "strand", "logFC_match", "FDR",
        "multiple_peaks", "region_pos", "n_regions_dm"
    )] %>%
        dplyr::mutate(
            n_regions_not_dm = NA,
            peak_id = paste0("PS_", .data$peak_id)
        )
    gaps.df <- gaps.df[, c(
        "peak_id", "seqnames", "start", "end", "width",
        "strand", "logFC_match", "FDR", "multiple_peaks", "region_pos",
        "n_regions_dm", "n_regions_not_dm"
    )]

    df <- rbind(peaks.df, gaps.df)
    df <- df %>%
        dplyr::filter(n_regions_dm > 2) %>%
        .[order(.$FDR), ] %>%
        dplyr::mutate(rank_p = seq_len(nrow(.)))

    df <- df[, c(
        "peak_id", "seqnames", "start", "end", "width", "strand",
        "rank_p", "logFC_match", "FDR", "multiple_peaks", # "region_pos",
        "n_regions_dm", "n_regions_not_dm"
    )]
    df
}







#' New: identify peaks from dm results
#'
#' @description
#' `aggregate_peaks()` aggregates upregulated regions as outputted from `edgeR_results()`. These peaks represent the region that the transcription factor of interest bound in.
#'
#' Peaks are ranked using the theory behind `csaw::getBestTest()`. For each peak, the region with the lowest adjusted p-value is used as the FDR for the peak, and the corresponding logFC is used as the logFC. Using the logFC in this way is useful, as the average logFC may skew the actual data.
#' Peaks with a gap between them of <= 150 bp are combined into 1 peak, accounting for many of the small regions having little data.
#'
#' @param dm_results A data.frame of differential testing results for each GATC region as outputted from `edgeR_results()`
#' @param gap_size The maximum gap length to include in a peak that separates two significantly enriched regions (peaks). Default is 150, based on an average sequencing of 75bp.
#' @return A `data.frame` of peaks. Columns are as follows: peak_id (Unique peak identifier, used internally - PS indicates a single peak, PM indicates the peak was combined), seqnames, start, end, width, strand, rank_p, logFC_match, FDR, multiple_peaks (number of peaks, NA if 1), n_regions_dm, n_regions_not_dm
#' @export
aggregate_peaks <- function(dm_results, gap_size = 150) {
    if (!is.data.frame(dm_results)) {
        stop("Must have data frame of differential_testing results from `edgeR_results")
    }
    df_aa <- dm_results
    df_a <- df_aa %>% dplyr::mutate(
        number = seq_len(dplyr::n()),
        trial = unsplit(lapply(split(.[, "meth_status"], .data$seqnames), function(x) {
            sequence(rle(x)$lengths)
        }), .data$seqnames),
        trial = ifelse(dplyr::lead(.data$trial) == .data$trial, 0, .data$trial),
        multiple = ifelse(.data$trial == 0, FALSE, TRUE)
    )
    df_1 <- peak_helper(df_a, TRUE, meth_status)
    df_2 <- peak_helper(df_a, FALSE, NULL)
    df_3 <- df_aa %>%
        dplyr::mutate(
            peak_id = df_1[match(df_aa$Position, df_1$Position), "seq"],
            not_consec_dm = df_2[match(df_aa$Position, df_2$Position), "seq"]
        ) %>%
        dplyr::group_by(.data$peak_id, .data$not_consec_dm) %>%
        dplyr::mutate(
            n_regions_dm = dplyr::n(),
            distance_dm = sum(.data$width),
            dm_start = ifelse(dplyr::row_number() == 1, .data$start, NA),
            dm_end = ifelse(dplyr::row_number() == dplyr::n(), .data$end, NA)
        ) %>%
        tidyr::fill(dm_start) %>%
        dplyr::ungroup() %>%
        tidyr::fill(dm_end, .direction = "up")

    peaks <- df_3 %>%
        dplyr::group_by(.data$peak_id) %>%
        dplyr::mutate(
            FDR = min(.data$adjust.p),
            logFC_match = ifelse(.data$adjust.p == .data$FDR, .data$logFC, NA)
        ) %>%
        tidyr::fill(logFC_match, .direction = "downup") %>%
        dplyr::ungroup() %>%
        dplyr::distinct(
            .data$seqnames, .data$dm_start, .data$dm_end, .data$meth_status,
            .data$peak_id, .data$n_regions_dm, .data$logFC_match, .data$FDR
        ) %>%
        dplyr::filter(
            !is.na(.data$peak_id),
            # .data$n_regions_dm != 2,
            .data$meth_status == "Upreg"
        ) %>%
        .[order(.$FDR), ] %>%
        dplyr::mutate(rank_p = seq_len(dplyr::n())) %>%
        .[order(.$peak_id), ] %>%
        stats::setNames(c(colnames(.[1]), "start", "end", colnames(.[4:9]))) %>%
        plyranges::as_granges() %>%
        data.frame() %>%
        dplyr::group_by(.data$seqnames) %>%
        dplyr::mutate(gap = dplyr::lead(.data$start) - .data$end) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(number = seq_len(nrow(.))) %>%
        data.frame()

    gaps <- gaps_fn_new(peaks, gap_size)
    if (nrow(gaps) == 0) {
        peaks_new <- peaks %>%
            data.frame() %>%
            dplyr::mutate(
                multiple_peaks = NA,
                n_regions_not_dm = 0
            )
        peaks_new <- order_peaks(peaks_new)
        return(peaks_new)
    }
    gaps <- gaps %>%
        dplyr::mutate(multiple_peaks = .data$n_peaks)
    peaks_no_gap <- peaks %>%
        dplyr::filter(!.data$peak_id %in% as.double(unlist(strsplit(gaps$id, ",")))) %>%
        dplyr::mutate(
            multiple_peaks = NA,
            n_regions_not_dm = 0
        )
    gaps_regions <- plyranges::find_overlaps_within(
        plyranges::as_granges(df_aa),
        plyranges::as_granges(gaps)
    ) %>%
        data.frame()
    gaps_regions <- gaps_regions %>%
        dplyr::mutate(is_dm = ifelse(.data$meth_status == "Upreg", TRUE, FALSE)) %>%
        dplyr::group_by(.data$Pos) %>%
        dplyr::mutate(
            n_regions = dplyr::n(),
            FDR = min(.data$adjust.p),
            logFC_match = ifelse(.data$adjust.p == .data$FDR, .data$logFC, NA)
        ) %>%
        tidyr::fill(logFC_match, .direction = "downup") %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$Pos, .data$is_dm) %>%
        dplyr::mutate(n_dm = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            n_regions_dm = ifelse(.data$is_dm == TRUE, .data$n_dm, NA),
            n_regions_not_dm = ifelse(.data$is_dm == FALSE, .data$n_dm, NA)
        ) %>%
        dplyr::group_by(.data$Pos) %>%
        tidyr::fill(n_regions_dm) %>%
        tidyr::fill(n_regions_not_dm, .direction = "downup") %>%
        dplyr::ungroup() %>%
        data.frame()

    gaps$logFC_match <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "logFC_match"]
    gaps$FDR <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "FDR"]
    gaps$n_regions <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions"]
    gaps$n_regions_dm <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions_dm"]
    gaps$n_regions_not_dm <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions_not_dm"]

    gaps <- gaps %>%
        dplyr::mutate(peak_id = .data$first) %>%
      .[, c("peak_id", "seqnames", "start", "end", "width", "strand",
            "logFC_match", "FDR", "multiple_peaks", "n_regions_dm", "n_regions_not_dm")]
    peaks_no_gap <- peaks_no_gap[, colnames(gaps)]

    peaks_new <- rbind(peaks_no_gap, gaps) %>%
        data.frame()
    peaks_new <- dplyr::filter(peaks_new, n_regions_dm != 2)
    peaks_new <- order_peaks(peaks_new)
    peaks_new
}



gaps_fn_new <- function(df, gap_size = 150) {
    gaps_work <- rbind(c(1, 1, 1, 1)) %>%
        data.frame() %>%
        stats::setNames(c("seqnames", "id", "start", "end")) %>%
        dplyr::mutate(seqnames = as.character(.data$seqnames))
    end <- c()
    gaps <- df
    gaps <- gaps %>%
        dplyr::mutate(
            gap = dplyr::coalesce(.data$gap, 100000),
            gap = ifelse(.data$gap < 0, 10000, .data$gap)
        ) %>%
        dplyr::group_by(.data$seqnames) %>%
        dplyr::mutate(total_peaks = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(gap = ifelse(.data$total_peaks == 1, -5, .data$gap)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$gap > 0) %>%
        .[, !colnames(.) == "total_peaks"]
    if (nrow(gaps) == 0) {
        return(gaps)
    }
    gaps_check <- gaps %>%
        dplyr::mutate(gap2 = ifelse(.data$gap <= 150 & .data$width <= 10000, TRUE, FALSE)) %>%
        dplyr::filter(.data$gap2 != FALSE) %>%
        .[, !colnames(.) == "gap2"]

    if (nrow(gaps_check) == 0) {
        return(gaps_check)
    }

    i <- 1
    number <- c(seq_len(nrow(gaps)))
    j <- i
    for (i in number) {
        if (gaps[i, ]$gap <= gap_size && gaps[i, ]$width < 10000) {
            start <- gaps[i, ]$start
            id <- gaps[i, ]$peak_id
            for (j in (i + 1):nrow(gaps)) {
                end <- c(end, gaps[j, ]$end)
                number <- c(number, gaps[j, ]$number)
                id <- paste0(id, ",", gaps[j, ]$peak_id)
                if (gaps[j, ]$gap > 150) {
                    break
                }
            }
            gaps_work[nrow(gaps_work) + 1, ] <- c(gaps[utils::tail(i, 1), ]$seqnames, id, start, utils::tail(end, 1))
            if (utils::tail(number, 1) == 100) {
                break
            }
        }
    }
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
        dplyr::mutate(n_peaks = stringr::str_count(id, ",") + 1)
    gaps_work <- gaps_work %>%
        dplyr::mutate(
            first = substr(id, 1, regexpr(",", id) - 1),
            first = as.double(.data$first)
        ) %>%
        data.frame()
    gaps_work$seqnames <- gaps[match(gaps_work$first, gaps$peak_id), "seqnames"]
    gaps_work %>%
        plyranges::as_granges() %>%
        data.frame() %>%
        dplyr::mutate(Pos = paste0(.data$seqnames, "-", .data$start))
}

peak_helper <- function(df_a, multiple, meth_status) {
    df <- df_a %>%
        dplyr::filter(multiple == {{ multiple }}) %>%
        dplyr::group_by(.data$seqnames, {{ meth_status }}) %>%
        dplyr::mutate(
            swap = ifelse(dplyr::lag(.data$number) != (.data$number - 1), 1, 0),
            swap = dplyr::coalesce(.data$swap, 1)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$swap) %>%
        dplyr::mutate(
            seq = ifelse(.data$swap == 1, seq_len(nrow(.)), 0),
            seq = ifelse(.data$seq == 0, NA, .data$seq)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by({{ meth_status }}) %>%
        tidyr::fill(seq) %>%
        dplyr::ungroup() %>%
        data.frame()
    df
}

order_peaks <- function(peaks) {
    peaks_new <- peaks %>%
        .[order(.$FDR), ] %>%
        dplyr::mutate(rank_p = seq_len(nrow(.))) %>%
      .[, c("peak_id", "seqnames", "start", "end", "width", "strand", "rank_p",
            "logFC_match", "FDR", "multiple_peaks", "n_regions_dm", "n_regions_not_dm")] %>%
        dplyr::mutate(peak_id = ifelse(is.na(.data$multiple_peaks), paste0("PS_", .data$peak_id), paste0("PM_", .data$peak_id)))
    row.names(peaks_new) <- NULL
    peaks_new
}
