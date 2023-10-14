
#' OLD:Identify peaks by aggregating differentially methylated regions
#'
#' `aggregate_peaks_old` aggregates differentially methylated GATC regions into peaks. These peaks represent the region that the gene of interest bound in.
#'
#' @param dm_results data.frame of differential methylation results obtained from [edgeR_results]
#'
#' @return data.frame of peaks.
#' Columns are
#' * seqnames,
#' * start,
#' * end,
#' * width,
#' * consec_dm (peak identifier),
#' * n_regions_dm (n of dm regions within the peak),
#' * ave_logFC (average logFC calculated across regions),
#' * ave_pVal (average adjusted P Value calculated across the regions),
#' * rank_p (based on average p value).
#' * gap (gap in bp to the next peak),
#' * number (peak number in order of position)
#' @export
#need to rename fn - aggregatePeaks
aggregate_peaks_old <- function(dm_results) {
  if(!is.data.frame(dm_results)) {
    stop("Must have data frame of differential_testing results from `edgeR_results")
  }
  results <- dm_results
  df_a <- results %>%
      dplyr::mutate(number = 1:nrow(.),
                    trial = unsplit(lapply(split(.[,"meth_status"], .$seqnames), function(x) {sequence(rle(x)$lengths)}), .$seqnames),
                    trial = ifelse(dplyr::lead(trial) == trial, 0, trial),
                    multiple = ifelse(trial == 0, FALSE, TRUE))
  df_1 <- df_a %>%
      dplyr::filter(multiple == TRUE) %>%
      dplyr::group_by(seqnames, meth_status) %>%
      dplyr::mutate(swap = ifelse(dplyr::lag(number) != (number - 1), 1, 0),
                    swap = dplyr::coalesce(swap, 1)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(swap) %>%
      dplyr::mutate(seq = ifelse(swap == 1, 1:nrow(.), 0),
                    seq = ifelse(seq == 0, NA, seq)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(meth_status) %>%
      tidyr::fill(seq) %>%
      dplyr::ungroup() %>%
      data.frame()
  df_2 <- df_a %>%
      dplyr::filter(multiple == FALSE) %>%
      dplyr::group_by(seqnames) %>%
      dplyr::mutate(swap = ifelse(dplyr::lag(number) != (number - 1), 1, 0),
                    swap = dplyr::coalesce(swap, 1)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(swap) %>%
      dplyr::mutate(seq = ifelse(swap == 1, 1:nrow(.), 0),
                    seq = ifelse(seq == 0, NA, seq)) %>%
      dplyr::ungroup() %>%
      tidyr::fill(seq) %>%
      dplyr::ungroup() %>%
      data.frame()
  df_3 <- results %>%
      dplyr::mutate(consec_dm = df_1[match(.$Position, df_1$Position), "seq"],
                    not_consec_dm = df_2[match(.$Position, df_2$Position), "seq"]) %>%
      dplyr::group_by(consec_dm, not_consec_dm) %>%
      dplyr::mutate(n_regions_dm = dplyr::n(),
                    distance_dm = sum(width),
                    dm_start = ifelse(dplyr::row_number() == 1, start, NA),
                    dm_end = ifelse(dplyr::row_number() == dplyr::n(), end, NA)) %>%
      tidyr::fill(dm_start) %>%
      dplyr::ungroup() %>%
      tidyr::fill(dm_end, .direction = "up")

  peaks <- df_3 %>%
      dplyr::group_by(consec_dm) %>%
      dplyr::mutate(ave_logFC = mean(logFC),
                    ave_pVal = mean(adjust.p)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(seqnames, dm_start, dm_end, meth_status, consec_dm, n_regions_dm, ave_logFC, ave_pVal) %>%
      dplyr::filter(!is.na(consec_dm), n_regions_dm != 2, meth_status == "Upreg") %>%
      .[order(.$ave_pVal),] %>%
      dplyr::mutate(rank_p = 1:nrow(.)) %>%
      .[order(.$consec_dm),] %>%
      stats::setNames(c(colnames(.[1]), "start", "end", colnames(.[4:9]))) %>%
      plyranges::as_granges() %>%
      data.frame() %>%
      dplyr::group_by(seqnames) %>%
      dplyr::mutate(gap = dplyr::lead(start) - end) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(number = 1:nrow(.)) %>%
      data.frame()
  peaks
}


#' New: identify peaks from de results
#'
#' @description
#' `aggregate_peaks()` aggregates upregulated regions as outputted from `edgeR_results()`. These peaks represent the region that the transcription factor of interest bound in.
#'
#' The average logFC and adjusted p-value is calculated for each peak, and peaks are ranked by their new p-value. Peaks with a gap between them of <= 150 bp are combined into 1 peak, accounting for many of the small regions having little data.
#'
#' @param dm_results A data.frame of differential testing results for each GATC region as outputted from `edgeR_results()`
#' @return A data.frame of peaks. Columns are as follows: peak_id (Unique peak identifier, used internally - PS indicates a single peak, PM indicates the peak was combined), seqnames, start, end, width, strand, rank_p, ave_logFC, ave-pVal, multiple_peaks (number of peaks, NA if 1), n_regions_dm, n_regions_not_dm
#' @export
#'
#' @examples
#' edgeR_results <- random_edgeR_result()
#'
#' aggregate_peaks(edgeR_results)
aggregate_peaks <- function(dm_results) {
  if(!is.data.frame(dm_results)) {
    stop("Must have data frame of differential_testing results from `edgeR_results")
  }
  df_aa <- dm_results
  df_a <- df_aa %>% dplyr::mutate(number = 1:dplyr::n(),
                                  trial = unsplit(lapply(split(.[,"meth_status"], .$seqnames), function(x) {sequence(rle(x)$lengths)}), .$seqnames),
                                  trial = ifelse(dplyr::lead(trial) == trial, 0, trial),
                                  multiple = ifelse(trial == 0, FALSE, TRUE))
  df_1 <- peak_helper(df_a, TRUE, meth_status)
  df_2 <- peak_helper(df_a, FALSE, NULL)
  df_3 <- df_aa %>%
    dplyr::mutate(peak_id = df_1[match(df_aa$Position, df_1$Position), "seq"],
                  not_consec_dm = df_2[match(df_aa$Position, df_2$Position), "seq"]) %>%
    dplyr::group_by(.data$peak_id, .data$not_consec_dm) %>%
    dplyr::mutate(n_regions_dm = dplyr::n(),
                  distance_dm = sum(.data$width),
                  dm_start = ifelse(dplyr::row_number() == 1, .data$start, NA),
                  dm_end = ifelse(dplyr::row_number() == dplyr::n(), .data$end, NA)) %>%
    tidyr::fill(.data$dm_start) %>%
    dplyr::ungroup() %>%
    tidyr::fill(.data$dm_end, .direction = "up")

  peaks <- df_3 %>%
    dplyr::group_by(.data$peak_id) %>%
    dplyr::mutate(ave_logFC = mean(.data$logFC),
                  ave_pVal = mean(.data$adjust.p)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$seqnames, .data$dm_start, .data$dm_end,
                    .data$meth_status, .data$peak_id, .data$n_regions_dm, .data$ave_logFC, .data$ave_pVal) %>%
    dplyr::filter(!is.na(.data$peak_id),
                  .data$n_regions_dm != 2,
                  .data$meth_status == "Upreg") %>%
    .[order(.$ave_pVal),] %>%
    dplyr::mutate(rank_p = 1:dplyr::n()) %>%
    .[order(.$peak_id),] %>%
    stats::setNames(c(colnames(.[1]), "start", "end", colnames(.[4:9]))) %>%
    plyranges::as_granges() %>%
    data.frame() %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::mutate(gap = dplyr::lead(.data$start) - .data$end) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(number = 1:nrow(.)) %>%
    data.frame()

  gaps <- gaps_fn_new(peaks)
  if(nrow(gaps) == 0) {
    peaks_new <- peaks %>%
      data.frame() %>%
      dplyr::mutate(multiple_peaks = NA,
                    n_regions_not_dm = 0) %>%
      .[,c(1:5,7,9:10,14,8,15)]
    peaks_new <- order_peaks(peaks_new)
    return(peaks_new)
  }
  gaps <- gaps %>%
    dplyr::mutate(multiple_peaks = .data$n_peaks)
  peaks_no_gap <- peaks %>%
    dplyr::filter(!peak_id %in% as.double(unlist(strsplit(gaps$id, ",")))) %>%
    dplyr::mutate(multiple_peaks = NA,
                  n_regions_not_dm = 0)
  gaps_regions <- plyranges::find_overlaps_within(plyranges::as_granges(df_aa),
                                                  plyranges::as_granges(gaps)) %>%
    data.frame()
  gaps_regions <- gaps_regions %>%
    dplyr::mutate(is_de = ifelse(.data$meth_status == "Upreg", TRUE, FALSE)) %>%
    dplyr::group_by(.data$Pos) %>%
    dplyr::mutate(n_regions = dplyr::n(),
                  ave_logFC = mean(.data$logFC),
                  ave_pVal = mean(.data$adjust.p)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Pos, .data$is_de) %>%
    dplyr::mutate(n_de = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n_regions_dm = ifelse(.data$is_de == TRUE, .data$n_de, NA),
                  n_regions_not_dm = ifelse(.data$is_de == FALSE, .data$n_de, NA)) %>%
    dplyr::group_by(.data$Pos) %>%
    tidyr::fill(.data$n_regions_dm) %>%
    tidyr::fill(.data$n_regions_not_dm, .direction = "downup") %>%
    dplyr::ungroup() %>%
    data.frame()

  gaps$ave_logFC <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "ave_logFC"]
  gaps$ave_pVal <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "ave_pVal"]
  gaps$n_regions <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions"]
  gaps$n_regions_dm <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions_dm"]
  gaps$n_regions_not_dm <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions_not_dm"]

  gaps <- gaps %>%
    dplyr::mutate(peak_id = .data$first) %>%
    .[,c(1:5,16,11:12,10,14,15)]
  peaks_no_gap <- peaks_no_gap[,c(1:5,7,9:10,14,8,15)]

  peaks_new <- rbind(peaks_no_gap, gaps) %>%
    data.frame()
  peaks_new <- order_peaks(peaks_new)
  peaks_new

}


#' Add de_results to full region df
#'
#' Used within aggregate_peaks to add in the regions that were excluded from edgeR analysis for low counts
#' Is also required for some plotting fns
#'
#' @param de_results as outputted from [edgeR_results()]
#' @param regions data.frame of regions, default is regions_gatc_drosophila_dm6
#'
#' @return full data.frame of regions with added information about the de results;
#' * de - 1,0,-1,NA ;
#' * logFC: 0 if de is NA ;
#' * adjust.p: 1 if de is NA :
#' * meth_status: Upreg, No_sig, Downreg, Not_included
#' @export
add_de <- function(de_results, regions=regions_gatc_drosophila_dm6) {
  if(!is.data.frame(de_results)) {
    stop("Must have data frame of differential_testing results from `edgeR_results")
  }
  if(!is.data.frame(regions)) {
    stop("Regions must be a data.frame")
  }
  results <- de_results
  df <- regions %>%
      dplyr::mutate(seqnames = paste0("chr", seqnames), number = 1:nrow(.))
  df$de <- results[match(df$Position, row.names(results)), "de"]
  df$logFC <- results[match(df$Position, row.names(results)), "logFC"]
  df$adjust.p <- results[match(df$Position, row.names(results)), "adjust.p"]
  df <- df %>%
      dplyr::mutate(meth_status = dplyr::case_when(is.na(.data$de) ~ "Not_included",
                                                   .data$de == 1 ~ "Upreg",
                                                   .data$de == -1 ~ "Downreg",
                                                   TRUE ~ "No_sig"))
  df <- df %>%
      dplyr::mutate(logFC = dplyr::coalesce(.data$logFC, 0), adjust.p = dplyr::coalesce(.data$adjust.p, 1))
  df
}



gaps_fn_new <- function(df) {
  gaps_work <- rbind(c(1, 1, 1, 1)) %>%
    data.frame() %>%
    stats::setNames(c("seqnames", "id", "start", "end")) %>%
    dplyr::mutate(seqnames = as.character(seqnames))
  end <- c()
  gaps <- df
  gaps <- gaps %>%
    dplyr::mutate(gap = dplyr::coalesce(.data$gap, 100000),
                  gap = ifelse(.data$gap < 0, 10000, .data$gap)) %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::mutate(total_peaks = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gap = ifelse(.data$total_peaks == 1, -5, .data$gap)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$gap > 0) %>%
    .[,!colnames(.) == "total_peaks"]
  if(nrow(gaps) == 0) {
    return(gaps)
  }
  gaps_check <- gaps %>%
    dplyr::mutate(gap2 = ifelse(.data$gap <= 150 & .data$width <= 10000, TRUE, FALSE)) %>%
    dplyr::filter(.data$gap2 != FALSE) %>%
    .[,!colnames(.) == "gap2"]

  if(nrow(gaps_check) == 0) {
    return(gaps_check)
  }

  i <- 1
  number <- c(1:nrow(gaps))
  j <- i
  for(i in number) {
    if(gaps[i,]$gap <= 150 && gaps[i,]$width < 10000) {
      start <- gaps[i,]$start
      id <- gaps[i,]$peak_id
      for(j in (i+1):nrow(gaps)) {
        end <- c(end, gaps[j,]$end)
        number <- c(number, gaps[j,]$number)
        id <- paste0(id, ",", gaps[j,]$peak_id)
        if(gaps[j,]$gap > 150) {
          break
        }
      }
      gaps_work[nrow(gaps_work) + 1,] <- c(gaps[tail(i, 1),]$seqnames, id, start, tail(end, 1))
      if(tail(number, 1) == 100) {
        break
      }
    }
  }
  gaps_work <- gaps_work[2:nrow(gaps_work),] %>%
    dplyr::group_by(.data$end) %>%
    dplyr::mutate(number = 1:dplyr::n(),
                  start = as.numeric(.data$start),
                  end = as.numeric(.data$end)) %>%
    dplyr::filter(.data$number == 1) %>%
    data.frame() %>%
    .[,colnames(.) != "number"] %>%
    dplyr::mutate(n_peaks = stringr::str_count(id, ",") + 1)
  gaps_work <- gaps_work %>%
    dplyr::mutate(first = substr(id, 1, regexpr(",", id) - 1),
                  first = as.double(.data$first)) %>%
    data.frame()
  gaps_work$seqnames <- gaps[match(gaps_work$first, gaps$peak_id), "seqnames"]
  gaps_work %>%
    plyranges::as_granges() %>%
    data.frame() %>%
    dplyr::mutate(Pos = paste0(seqnames, "-", start))

}

peak_helper <- function(df_a, multiple, meth_status) {
  df <- df_a %>%
    dplyr::filter(multiple == {{multiple}}) %>%
    dplyr::group_by(seqnames, {{meth_status}}) %>%
    dplyr::mutate(swap = ifelse(dplyr::lag(.data$number) != (.data$number - 1), 1, 0),
                  swap = dplyr::coalesce(.data$swap, 1)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$swap) %>%
    dplyr::mutate(seq = ifelse(.data$swap == 1, 1:nrow(.), 0),
                  seq = ifelse(.data$seq == 0, NA, seq)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by({{meth_status}}) %>%
    tidyr::fill(.data$seq) %>%
    dplyr::ungroup() %>%
    data.frame()
  df
}

order_peaks <- function(peaks) {
  peaks_new <- peaks %>%
    .[order(.$ave_pVal),] %>%
    dplyr::mutate(rank_p = 1:nrow(.)) %>%
    .[,c(6,1:5,12,7:11)] %>%
    .[order(.$peak_id),] %>%
    dplyr::mutate(peak_id = ifelse(is.na(.data$multiple_peaks), paste0("PS_", peak_id), paste0("PM_", peak_id)))
  row.names(peaks_new) <- NULL
  peaks_new
}
