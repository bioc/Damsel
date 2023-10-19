


#' Create example regions
#'
#' @param size number of rows to create
#'
#' @return example data.frame with output similar to `gatc_track()$regions`
#' @export
random_regions <- function(size=50) {
  df <- list(start = 50, end = 85, width = 36) %>% data.frame()
  size_n <- size - 1
  random_width <- floor(stats::runif(size_n, 5, 1000))
  for(i in 1:size_n) {
    new_start <- df[i, "end"] + 1
    df[nrow(df) + 1,] <- list(start = new_start, end = (new_start + random_width[i] - 1), width = random_width[i])
  }
  df$seqnames <- "2L"
  df$strand <- "*"
  row.names(df) <- NULL
  df$Position <- paste0("chr", df$seqnames, "-", df$start)
  df[,c("Position", "seqnames", "start", "end", "width", "strand")]
}



#' Create example counts
#'
#' @param size number of rows to create
#'
#' @return example data.frame of counts similar to `process_bams()`
#' @export
random_counts <- function(size=50) {
  counts <- random_regions(size)
  size <- nrow(counts)
  counts$Dam_1.bam <- stats::rnorm(size, 100)
  counts$Fusion_1.bam <- stats::rnorm(size, 400)
  counts$Dam_2.bam <- counts$Dam_1.bam + 7
  counts$Fusion_2.bam <- counts$Fusion_1.bam - 2
  counts$seqnames <- paste0("chr", counts$seqnames)
  counts
}


#' Create example edgeR results
#'
#' @param size number of rows to create
#'
#' @return example data.frame of edgeR results, output similar to `edgeR_results()`
#' @export
random_edgeR_results <- function(size=50) {
  results <- random_regions(size)
  results$seqnames <- paste0("chr", results$seqnames)
  results$number <- 1:size
  dm_options <- c(0, 1, NA, -1)
  peak <- rep(1, each = 4)
  zero <- rep(0, each = 2)
  pairs <- rep(dm_options, each = 2)
  threes <- rep(dm_options, each = 3)
  all <- c(dm_options, dm_options, dm_options, dm_options, peak, pairs, peak, zero, threes, peak)
  results$dm <- all
  results$logFC <- dplyr::case_when(results$dm == 1 ~ runif(1, 1, 5),
                                    results$dm == 0 ~ runif(1, -1.5, 1.5),
                                    results$dm == -1 ~ runif(1, -7, -1),
                                    TRUE ~ 0)
  results$adjust.p <- dplyr::case_when(results$dm == 1 ~ runif(1, 5.58e-07, 0.05),
                                       results$dm == 0 ~ runif(1, 9.4e-02, 1),
                                       results$dm == -1 ~ runif(1, 2.93e-06, 0.05),
                                       TRUE ~ 1)
  results$meth_status <- dplyr::case_when(results$dm == 1 ~ "Upreg",
                                          results$dm == 0 ~ "No_signal",
                                          results$dm == -1 ~ "Downreg",
                                          TRUE ~ "Not_included")
  results
}
