##function
find_peaks <- function(de_results, regions = regions) {
  df <- regions %>% dplyr::mutate(sig = de_results[match(.$Position, row.names(de_results)), "de"],
                           logFC = de_results[match(.$Position, row.names(de_results)), "logFC"],
                           pVal = de_results[match(.$Position, row.names(de_results)), "PValue"])
  df <- df %>% dplyr::mutate(sig = dplyr::coalesce(sig, 0), logFC = dplyr::coalesce(logFC, 0), pVal = dplyr::coalesce(pVal, 1))
  df$tag <- df$start

  start <- -1
  end <- -1
  logFC <- 0
  size <- 0
  sig <- 0

  lapseOneTag <- FALSE
  peaks <- rbind(1:7) %>% data.frame()
  for (i in seq(1, nrow(df), by = 2)) {
    if (i + 3 < nrow(df) && df[i + 3,]$pVal < 0.05 && df[i,]$pVal < 0.05 &&
        df[i + 3,]$logFC > 0 && df[i,]$logFC > 0) {
      if (df[i,]$sig > 0) {
        sig <- sig + 1
      }
      if (df[i + 3,]$sig > 0) {
        sig <- sig + 1
      }
      size <- size + 2
      logFC <- logFC + df[i,]$logFC + df[i + 3,]$logFC
      if (start < 0) {
        start <- i
      }
      end <- i + 3
      lapseOneTag <- TRUE
    } else if (lapseOneTag) {
      lapseOneTag <- FALSE
    } else if (start >= 0) {
      s <-  df[start,]$tag
      e <- df[end,]$tag
      m <- (s + e) / 2
      unmeth <- 0
      for (j in seq(end + 3, nrow(df), by = 1)) {
        x <- df[j,]$tag
        if (x - m > 500) {
          break
        }
        if (!df[j,]$sig && !df[j - 1,]$sig) {
          unmeth <- unmeth + 1
        }
      }
      if(start == 1) {
        for (j in seq(1, start - 1, by= -1)) {
          x <- df[j,]$tag
          if (m - x > 500) {
            break
          }
          if (!df[j,]$sig && !df[j + 1,]$sig) {
            unmeth <- unmeth + 1
          }
        } } else {
          for (j in seq(start - 1, 1, by= -1)) {
            x <- df[j,]$tag
            if (m - x > 500) {
              break
            }
            if (!df[j,]$sig && !df[j + 1,]$sig) {
              unmeth <- unmeth + 1
            }
          }
        }

      ll <- logFC / size
      output <- c(df[start,]$tag, s, e, size, (size / (unmeth + size)), ll, sig)
      peaks[nrow(peaks) + 1, ] <- output

      last <- 0
      skip <- 0
      logFC <- 0
      start <- -1
      end <- -1
      size <- 0
      sig <- 0

    }
  }

  peaks <- peaks[2:nrow(peaks), 2:ncol(peaks)]
  colnames(peaks) <- c("start", "end", "n_regions", "any_unmeth", "aveLogFC", "n_de_regions")
  peaks$big = peaks$n_regions> 2 | peaks$any_unmeth==1

  peaks$number <- 1:nrow(peaks)

  seqnames_pos <- peaks %>%
    dplyr::mutate(prev_ = lag(start), start_diff = start-prev_) %>%
    dplyr::filter(start_diff < 0) %>% .[, "number"]
  seqnames_pos <- c(1, seqnames_pos)
  seqnames_pos <- c(seqnames_pos, (nrow(peaks) + 1))
  seqnames_distinct <- df$seqnames %>% unique()

  seqnames_vec <- character()
  for(i in 1:length(seqnames_distinct)) {
    out <- replicate(n = (seqnames_pos[i + 1] - seqnames_pos[i]), seqnames_distinct[i])
    seqnames_vec <- c(seqnames_vec, out)
  }

  peaks$seqnames <- seqnames_vec
  peaks
}

