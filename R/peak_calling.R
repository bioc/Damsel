
#' Identify peaks by aggregating differentially methylated regions
#'
#' `aggregate_peaks` aggregates differentially methylated GATC regions into peaks. These peaks represent the region that the gene of interest bound in.
#'
#' @param dm_results data frame of differential methylation results obtained from [edgeR_results]
#' @param regions data frame of GATC regions. Default is GATC regions from Drosophila melanogaster - dm6.
#'
#' @return data frame of peaks.
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
#'
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' dge <- edgeR_set_up(counts.df)
#' de_results <- edgeR_results(dge, p.value = 0.05, lfc = 1)
#'
#' peaks <- aggregate_peaks(de_results, regions = regions_gatc_drosophila_dm6)
#' head(peaks)
#need to rename fn - aggregatePeaks
aggregate_peaks <- function(dm_results, regions=regions_gatc_drosophila_dm6) {
  if(!is.data.frame(dm_results)) {
    stop("Must have data frame of differential_testing results from `edgeR_results")
  }
  if(!is.data.frame(regions)) {
    stop("Regions must be a data.frame")
  }
  if(missing(regions)) {
    message("Default of drosophila dm6 regions used")
  }
  results <- add_de(de_results=dm_results, regions=regions)
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
#' @param dm_results data.frame as outputted from `edgeR_results()`
#' @param regions regions df. default is the provided data -
#' `regions_gatc_drosophila_dm6`
#'
#' @return data.frame of peaks
#' @export
#'
#' @examples
aggregate_peaks_new <- function(dm_results, regions=regions_gatc_drosophila_dm6) {
  if(!is.data.frame(dm_results)) {
    stop("Must have data frame of differential_testing results from `edgeR_results")
  }
  if(!is.data.frame(regions)) {
    stop("Regions must be a data.frame")
  }
  if(missing(regions)) {
    message("Default of drosophila dm6 regions used")
  }
  df_aa <- add_de(de_results=dm_results, regions=regions)
  df_a <- df_aa %>% dplyr::mutate(number = 1:dplyr::n(),
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
  df_3 <- df_aa %>%
    dplyr::mutate(peak_id = df_1[match(df_aa$Position, df_1$Position), "seq"],
                  not_consec_dm = df_2[match(df_aa$Position, df_2$Position), "seq"]) %>%
    dplyr::group_by(peak_id, not_consec_dm) %>%
    dplyr::mutate(n_regions_dm = dplyr::n(),
                  distance_dm = sum(width),
                  dm_start = ifelse(dplyr::row_number() == 1, start, NA),
                  dm_end = ifelse(dplyr::row_number() == dplyr::n(), end, NA)) %>%
    tidyr::fill(dm_start) %>%
    dplyr::ungroup() %>%
    tidyr::fill(dm_end, .direction = "up")

  peaks <- df_3 %>%
    dplyr::group_by(peak_id) %>%
    dplyr::mutate(ave_logFC = mean(logFC),
                  ave_pVal = mean(adjust.p)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(seqnames, dm_start, dm_end,
                    meth_status, peak_id, n_regions_dm, ave_logFC, ave_pVal) %>%
    dplyr::filter(!is.na(peak_id),
                  n_regions_dm != 2,
                  meth_status == "Upreg") %>%
    .[order(.$ave_pVal),] %>%
    dplyr::mutate(rank_p = 1:dplyr::n()) %>%
    .[order(.$peak_id),] %>%
    stats::setNames(c(colnames(.[1]), "start", "end", colnames(.[4:9]))) %>%
    plyranges::as_granges() %>%
    data.frame() %>%
    dplyr::group_by(seqnames) %>%
    dplyr::mutate(gap = dplyr::lead(start) - end) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(number = 1:nrow(.)) %>%
    data.frame()

  gaps <- gaps_fn_new(peaks)
  if(nrow(gaps) == 0) {
    peaks_new <- peaks %>%
      data.frame() %>%
      dplyr::mutate(multiple_peaks = NA,
                    n_regions_not_dm = 0) %>%
      .[,c(1:5,7,9:10,14,8,15)] %>%
      .[order(.$ave_pVal),] %>%
      dplyr::mutate(rank_p = 1:nrow(.)) %>%
      .[,c(6,1:5,12,7:11)] %>%
      .[order(.$peak_id),] %>%
      dplyr::mutate(peak_id = ifelse(is.na(multiple_peaks),
                                     paste0("PS_", peak_id),
                                     paste0("PM_", peak_id)))
    row.names(peaks_new) <- NULL
    return(peaks_new)
  }
  gaps <- gaps %>%
    dplyr::mutate(multiple_peaks = n_peaks)
  peaks_no_gap <- peaks %>%
    dplyr::filter(!peak_id %in% as.double(unlist(strsplit(gaps$id, ",")))) %>%
    dplyr::mutate(multiple_peaks = NA,
                  n_regions_not_dm = 0)
  gaps_regions <- plyranges::find_overlaps_within(plyranges::as_granges(df_aa),
                                                  plyranges::as_granges(gaps)) %>%
    data.frame()
  gaps_regions <- gaps_regions %>%
    dplyr::mutate(is_de = ifelse(meth_status == "Upreg", TRUE, FALSE)) %>%
    dplyr::group_by(Pos) %>%
    dplyr::mutate(n_regions = dplyr::n(),
                  ave_logFC = mean(logFC),
                  ave_pVal = mean(adjust.p)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Pos, is_de) %>%
    dplyr::mutate(n_de = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n_regions_dm = ifelse(is_de == TRUE, n_de, NA),
                  n_regions_not_dm = ifelse(is_de == FALSE, n_de, NA)) %>%
    dplyr::group_by(Pos) %>%
    tidyr::fill(n_regions_dm) %>%
    tidyr::fill(n_regions_not_dm, .direction = "downup") %>%
    dplyr::ungroup() %>%
    data.frame()

  gaps$ave_logFC <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "ave_logFC"]
  gaps$ave_pVal <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "ave_pVal"]
  gaps$n_regions <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions"]
  gaps$n_regions_dm <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions_dm"]
  gaps$n_regions_not_dm <- gaps_regions[match(gaps$Pos, gaps_regions$Pos), "n_regions_not_dm"]

  gaps <- gaps %>%
    dplyr::mutate(peak_id = first) %>%
    .[,c(1:5,16,11:12,10,14,15)]
  peaks_no_gap <- peaks_no_gap[,c(1:5,7,9:10,14,8,15)]

  peaks_new <- rbind(peaks_no_gap, gaps) %>%
    data.frame() %>% .[order(.$ave_pVal),] %>%
    dplyr::mutate(rank_p = 1:nrow(.)) %>%
    .[,c(6,1:5,12,7:11)] %>%
    .[order(.$peak_id),] %>%
    dplyr::mutate(peak_id = ifelse(is.na(multiple_peaks), paste0("PS_", peak_id), paste0("PM_", peak_id)))
  row.names(peaks_new) <- NULL
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
#'
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' dge <- edgeR_set_up(counts.df)
#' de_results <- edgeR_results(dge, p.value = 0.05, lfc = 1)
#'
#' de_results <- add_de(de_results, regions = regions_gatc_drosophila_dm6)
#' head(de_results)
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
      dplyr::mutate(meth_status = dplyr::case_when(is.na(de) ~ "Not_included",
                                                   de == 1 ~ "Upreg",
                                                   de == -1 ~ "Downreg",
                                                   TRUE ~ "No_sig"))
  df <- df %>%
      dplyr::mutate(logFC = dplyr::coalesce(logFC, 0), adjust.p = dplyr::coalesce(adjust.p, 1))
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
    dplyr::mutate(gap = dplyr::coalesce(gap, 100000),
                  gap = ifelse(gap < 0, 10000, gap)) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::mutate(gap = ifelse(dplyr::n() == 1, -5, gap)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(gap > 0)
  if(nrow(gaps) == 0) {
    return(gaps)
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
    dplyr::group_by(end) %>%
    dplyr::mutate(number = 1:dplyr::n(),
                  start = as.numeric(start),
                  end = as.numeric(end)) %>%
    dplyr::filter(number == 1) %>%
    data.frame() %>%
    .[,colnames(.) != "number"] %>%
    dplyr::mutate(n_peaks = stringr::str_count(id, ",") + 1)
  gaps_work <- gaps_work %>%
    dplyr::mutate(first = substr(id, 1, regexpr(",", id) - 1),
                  first = as.double(first)) %>%
    data.frame()
  gaps_work$seqnames <- gaps[match(gaps_work$first, gaps$peak_id), "seqnames"]
  gaps_work %>%
    plyranges::as_granges() %>%
    data.frame() %>%
    dplyr::mutate(Pos = paste0(seqnames, "-", start))

}
