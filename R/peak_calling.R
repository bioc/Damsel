
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
aggregate_peaks <- function(dm_results) {
  df_a <- dm_results %>%
            dplyr::mutate(number = 1:n(),
            trial = unsplit(lapply(split(.[,"de"], .$seqnames), function(x) {sequence(rle(x)$lengths)}), .$seqnames),
            trial = ifelse(lead(trial) == trial, 0, trial),
            multiple = ifelse(trial == 0, FALSE, TRUE))
  df_1 <- df_a %>%
            dplyr::filter(multiple == TRUE) %>%
            dplyr::group_by(seqnames, de) %>%
            dplyr::mutate(swap = ifelse(lag(number) != (number - 1), 1, 0),
                          swap = dplyr::coalesce(swap, 1)) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(swap) %>%
            dplyr::mutate(seq = ifelse(swap == 1, 1:nrow(.), 0),
                          seq = ifelse(seq == 0, NA, seq)) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(sig) %>%
            tidyr::fill(seq) %>%
            dplyr::ungroup() %>%
            data.frame()
  df_2 <- df_a %>%
            dplyr::filter(multiple == FALSE)  %>%
            dplyr::group_by(seqnames) %>%
            dplyr::mutate(swap = ifelse(lag(number) != (number - 1), 1, 0),
                          swap = dplyr::coalesce(swap, 1)) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(swap) %>%
            dplyr::mutate(seq = ifelse(swap == 1, 1:nrow(.), 0),
                          seq = ifelse(seq == 0, NA, seq)) %>%
            dplyr::ungroup() %>%
            tidyr::fill(seq) %>%
            dplyr::ungroup() %>%
            data.frame()
  df_3 <- dm_results %>%
            dplyr::mutate(consec_dm = df_1[match(df$Position, df_1$Position), "seq"],
                          not_consec_dm = df_2[match(df$Position, df_2$Position), "seq"]) %>%
            dplyr::group_by(consec_dm, not_consec_dm) %>%
            dplyr::mutate(n_regions_dm = n(),
                          distance_dm = sum(width),
                          dm_start = ifelse(row_number() == 1, start, NA),
                          dm_end = ifelse(row_number() == n(), end, NA)) %>%
            tidyr::fill(dm_start) %>%
            dplyr::ungroup() %>%
            tidyr::fill(dm_end, .direction = "up")

  peaks <- df_3 %>%
             dplyr::group_by(consec_dm) %>%
             dplyr::mutate(ave_logFC = mean(logFC),
                           ave_pVal = mean(pVal)) %>%
             dplyr::ungroup() %>%
             dplyr::distinct(seqnames, dm_start, dm_end, sig, consec_dm, n_regions_dm, ave_logFC, ave_pVal) %>%
             dplyr::filter(!is.na(consec_dm), n_regions_dm != 2, sig == 1) %>%
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
             dplyr::mutate(info = case_when(width <= 2000 ~ "significant binding",
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
