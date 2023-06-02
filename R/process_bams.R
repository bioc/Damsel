process_bams <- function(path_to_bams, regions) {
  #list of all bam files
  files <- list.files(path_to_bams, pattern = ".bam")
  #remove any bai files
  files <- files %>% as.data.frame() %>%
    stats::setNames(., "file") %>%
    mutate(bai = stringr::str_detect(.$file, "bai")) %>%
    filter(bai == FALSE)
  #add / to get path for each file
  files <- paste0(path_to_bams, "/", files$file)
  df <- regions
  #test format of seqnames and separate them
  seqnames_vec <- character()
  for(i in 1:length(files)) {
    is_error <- ifelse(sum(exomeCopy::countBamInGRanges(files[[i]], plyranges::as_granges(df))) == 0, "error", "fine")
    seqnames_vec <- c(seqnames_vec, is_error)
  }
  files <- as.data.frame(files) %>% setNames(., "file")
  files <- cbind(files, seqnames_vec)
  files_fine <- filter(files, seqnames_vec == "fine")$file
  files_error <- filter(files, seqnames_vec == "error")$file
  #seqnames format 2L etc
  for(i in 1:length(files_fine)) {
    if(length(files_fine) == 0) {
      break
    } else {
      raw_counts <- exomeCopy::countBamInGRanges(files_fine[[i]], plyranges::as_granges(df))
      df[, ncol(df) + 1] <- raw_counts
      colnames(df)[ncol(df)] <- paste0("raw_counts", files_fine[[i]])
    }
  }
  #seqnames format chr2L etc
  df <- df %>% mutate(seqnames = paste0("chr", seqnames))
  for(i in 1:length(files_error)) {
    if(length(files_error) == 0) {
      break
    } else {
      raw_counts <- exomeCopy::countBamInGRanges(files_error[[i]], plyranges::as_granges(df))
      df[, ncol(df) + 1] <- raw_counts
      colnames(df)[ncol(df)] <- paste0("raw_counts", files_error[[i]])
    }
  }

  df

}
