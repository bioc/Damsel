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

corr_heatmap <- function(df, method = "spearman") {
  corr_res <- cor(df[,7:ncol(df)], method = method)
  median_corr <- round(median(corr_res), 1)
  min_corr <- floor(min(corr_res)*10)/10
  corr_res <- round(corr_res, 2)
  # Use correlation between variables as distance and reorder
  dd <- as.dist((1-corr_res)/2)
  hc <- hclust(dd)
  corr_res <- corr_res[hc$order, hc$order]
  # upper triangle
  corr_res[lower.tri(corr_res)]<- NA
  # Melt the correlation matrix
  corr_res <- reshape2::melt(corr_res, na.rm = TRUE)
  #plot heatmap
  ggplot(corr_res, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red",
                         midpoint = median_corr, limit = c(min_corr,1), space = "Lab",
                         name = paste(stringr::str_to_title(method), "\nCorrelation", sep = " ")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
    coord_fixed()
}

corr_scatter <- function(df, sample_1, sample_2, method = "spearman") {
  ggpubr::ggscatter(df, x = sample_1, y = sample_2,
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = method)
}
