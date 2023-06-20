#main function
consec2_fn <- function(df) {
  df_1 <- df %>% mutate(trial = unsplit(lapply(split(.[,"meth_status"], .$seqnames), function(x) {sequence(rle(x)$lengths)}), .$seqnames)) %>%
    mutate(trial = ifelse(lead(trial) == trial, 0, trial)) %>%
    mutate(multiple = ifelse(trial == 0, F, T)) %>% filter(multiple == TRUE) %>%
    group_by(seqnames, meth_status) %>%
    mutate(swap = ifelse(lag(number) != (number - 1), 1, 0), swap = coalesce(swap, 1)) %>% ungroup() %>% group_by(swap) %>% mutate(seq = ifelse(swap == 1, 1:nrow(.), 0)) %>% mutate(seq = ifelse(seq == 0, NA, seq)) %>% ungroup() %>% group_by(meth_status) %>% fill(seq) %>% ungroup() %>% as.data.frame()
  df_2 <- df %>% mutate(trial = unsplit(lapply(split(.[,"meth_status"], .$seqnames), function(x) {sequence(rle(x)$lengths)}), .$seqnames)) %>%
    mutate(trial = ifelse(lead(trial) == trial, 0, trial)) %>%
    mutate(multiple = ifelse(trial == 0, F, T)) %>% filter(multiple == FALSE)  %>%
    group_by(seqnames) %>%
    mutate(swap = ifelse(lag(number) != (number - 1), 1, 0), swap = coalesce(swap, 1)) %>% ungroup() %>% group_by(swap) %>% mutate(seq = ifelse(swap == 1, 1:nrow(.), 0)) %>% mutate(seq = ifelse(seq == 0, NA, seq)) %>% ungroup() %>% fill(seq) %>% ungroup() %>% as.data.frame()
  df_3 <- df %>% mutate(consec_dm = df_1[match(df$Position, df_1$Position), "seq"], not_consec_dm = df_2[match(df$Position, df_2$Position), "seq"]) %>% group_by(consec_dm, not_consec_dm) %>% mutate(n_regions_dm = n(), distance_dm = sum(width), dm_start = ifelse(row_number() == 1, start, NA), dm_end = ifelse(row_number() == n(), end, NA)) %>% fill(dm_start) %>% ungroup() %>% fill(dm_end, .direction = "up")

  df <- df %>% mutate(size = case_when(width <= 500 ~ "Close", width > 500 & width <= 2000 ~ "Medium", TRUE ~ "Far"))

  df_4 <- df %>% mutate(trial = unsplit(lapply(split(.[,"size"], .$seqnames), function(x) {sequence(rle(x)$lengths)}), .$seqnames)) %>%
    mutate(trial = ifelse(lead(trial) == trial, 0, trial)) %>%
    mutate(multiple = ifelse(trial == 0, F, T)) %>% filter(multiple == TRUE) %>%
    group_by(seqnames, size) %>%
    mutate(swap = ifelse(lag(number) != (number - 1), 1, 0), swap = coalesce(swap, 1)) %>% ungroup() %>% group_by(swap) %>% mutate(seq = ifelse(swap == 1, 1:nrow(.), 0)) %>% mutate(seq = ifelse(seq == 0, NA, seq)) %>% ungroup() %>% group_by(size) %>% fill(seq) %>% ungroup() %>% as.data.frame()
  df_5 <- df %>% mutate(trial = unsplit(lapply(split(.[,"size"], .$seqnames), function(x) {sequence(rle(x)$lengths)}), .$seqnames)) %>%
    mutate(trial = ifelse(lead(trial) == trial, 0, trial)) %>%
    mutate(multiple = ifelse(trial == 0, F, T)) %>% filter(multiple == FALSE)  %>%
    group_by(seqnames) %>%
    mutate(swap = ifelse(lag(number) != (number - 1), 1, 0), swap = coalesce(swap, 1)) %>% ungroup() %>% group_by(swap) %>% mutate(seq = ifelse(swap == 1, 1:nrow(.), 0)) %>% mutate(seq = ifelse(seq == 0, NA, seq)) %>% ungroup() %>% fill(seq) %>% ungroup() %>% as.data.frame()
  df_3 <- df_3 %>% mutate(size = case_when(width <= 500 ~ "Close", width > 500 & width <= 2000 ~ "Medium", TRUE ~ "Far"), consec_size = df_4[match(df_3$Position, df_4$Position), "seq"], not_consec_size = df_5[match(df_3$Position, df_5$Position), "seq"]) %>% group_by(consec_size, not_consec_size) %>% mutate(n_regions_size = n(), distance_size = sum(width), type = ifelse(!is.na(consec_size), size, "Mix"), group_start = ifelse(row_number() == 1, start, NA), group_end = ifelse(row_number() == n(), end, NA)) %>% fill(group_start) %>% ungroup() %>% fill(group_end, .direction = "up")
  df_3
}
#run this beforehand
add_de <- function(regions, de_results) {
  df <- regions %>% mutate(number = 1:nrow(.))
  df$diff_meth <- de_results[match(df$Position, row.names(de_results)), "de"]
  df <- df %>% mutate(meth_status = case_when(is.na(diff_meth) ~ "Not_included", diff_meth == 1 ~ "Upreg", diff_meth == -1 ~ "Downreg", TRUE ~ "No_sig"))
  df
}


#jan's peak caller

peaks_jan_fn <- function(regions = regions_between_gatc_dm6) {
  system('python3 call_peaks.py keep lrt_sd.txt > peaks.txt')
  peaks <- read.table("peaks.txt")

  df <- peaks
  names(df) <- c('chr', 'start', 'end', "tags", 'pen', 'aveLogFC', 'sig')
  df$big <- df$tags> 2 | df$pen==1
  df$number <- 1:nrow(df)

  seqnames_pos <- df %>%
    mutate(prev_ = lag(start), start_diff = start-prev_) %>%
    filter(start_diff < 0) %>% .[, "number"]
  seqnames_pos <- c(1, seqnames_pos)
  seqnames_pos <- c(seqnames_pos, (nrow(df) + 1))
  seqnames_distinct <- regions$seqnames %>% unique()

  seqnames_vec <- character()
  for(i in 1:length(seqnames_pos)) {
    out <- replicate(n = (seqnames_pos[i + 1] - seqnames_pos[i]), seqnames_distinct[i])
    seqnames_vec <- c(seqnames_vec, out)
  }

  df$seqnames <- seqnames_vec
  df
}
