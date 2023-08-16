
#' Identify peaks by aggregating differentially methylated regions
#'
#' `aggregate_peaks` aggregates differentially methylated GATC regions into peaks. These peaks represent the region that the gene of interest bound in.
#'
#' @param dm_results data frame of differential methylation results obtained from [edgeR_results]
#'
#' @return data frame of peaks. columns are chromosome name, start, end, and width of peak, peak identifier, n of dm regions within the peak, gap in bp to the next peak, and peak rank (based on p value).
#' @export
#'
#' @examples
#need to rename fn - aggregatePeaks
aggregate_peaks <- function(dm_results, regions=regions_gatc_drosophila_dm6) {
  if(missing(dm_results) | !is.data.frame(dm_results)) {
    stop("Must have data frame of differential_testing results from `edgeR_results")
  }
  if(!is.data.frame(regions)) {
    stop("Regions must be a data.frame")
  }
  results <- add_de(de_results=dm_results, regions=regions_gatc_drosophila_dm6)
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
      dplyr::distinct(seqnames, dm_start, dm_end, meth_status, de, consec_dm, n_regions_dm, ave_logFC, ave_pVal) %>%
      dplyr::filter(!is.na(consec_dm), n_regions_dm != 2, de == 1) %>%
      .[order(.$ave_pVal),] %>%
      dplyr::mutate(rank_p = 1:nrow(.)) %>%
      .[order(.$consec_dm),] %>%
      stats::setNames(c(colnames(.[1]), "start", "end", colnames(.[4:6]))) %>%
      plyranges::as_granges() %>%
      data.frame() %>%
      dplyr::group_by(seqnames) %>%
      dplyr::mutate(gap = start - lag(end)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(number = 1:nrow(.))
  peaks <- peaks %>%
      dplyr::mutate(info = dplyr::case_when(width <= 2000 ~ "significant binding",
                                     width > 2000 && width <= 10000 ~ "less significant",
                                     width > 10000 ~ "unexpected width")) #some kind of message about peaks that are too big or something
#change 1:nrow(peaks) to i in seq_len or something
  for(i in 1:nrow(peaks)) {
    if(peaks[i,]$gap == 5) {
      peaks[i-1,] <- peaks[i-1,] %>% dplyr::mutate(end = peaks[i,]$end,
                                                   n_regions_dm = n_regions_dm + peaks[i,]$n_regions_dm,
                                                   width = width + peaks[i,]$width + 4)
      peaks <- peaks[-i,]
    }
  }
  peaks
}

add_de <- function(de_results, regions=regions_gatc_drosophila_dm6) {
  if(missing(de_results) | !is.data.frame(de_results)) {
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
      dplyr::mutate(meth_status = dplyr::case_when(is.na(de) ~ "Not_included",
                                                   de == 1 ~ "Upreg",
                                                   de == -1 ~ "Downreg",
                                                   TRUE ~ "No_sig"))
  df <- df %>%
      dplyr::mutate(logFC = dplyr::coalesce(logFC, 0), adjust.p = dplyr::coalesce(adjust.p, 1))
  df
}
