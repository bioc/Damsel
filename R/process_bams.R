process_bams <- function(path_to_bams, regions, cores=detectCores()) {
  #list of all bam files
  files <- list.files(path_to_bams, pattern = ".bam")
  #remove any bai files
  files <- files %>% data.frame() %>%
    stats::setNames(., "file") %>%
    dplyr::mutate(bai = stringr::str_detect(.$file, "bai")) %>%
    dplyr::filter(bai == FALSE)
  #add / to get path for each file
  path_to_bams <- ifelse(substring(path_to_bams, first = nchar(path_to_bams)) == "/", path_to_bams, paste0(path_to_bams, "/"))
  files <- paste0(path_to_bams, files$file)
  df <- regions
  #test format of seqnames and separate them
  one_count <- parallel::mclapply(files, function(x) {exomeCopy::countBamInGRanges(x, plyranges::as_granges(df))}, mc.cores = cores) %>% data.frame()
  colnames(one_count) <- files

  for(i in 1:ncol(one_count)) {
    if(sum(one_count[,i]) == 0) {
      one_count[,i] <- exomeCopy::countBamInGRanges(files[i], plyranges::as_granges(dplyr::mutate(df, seqnames = paste0("chr", seqnames))))
    }
  }

  df2 <- cbind(df, one_count)
  df2$seqnames <- paste0("chr", df2$seqnames)
  colnames(df2) <- sub(path_to_bams, "", colnames(df2))
  df2

}


corr_heatmap <- function(df, method = "spearman") {
  corr_res <- cor(df[,grepl("bam", colnames(df))], method = method)
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
  ggplot::ggplot(corr_res, aes(Var2, Var1, fill = value)) +
    ggplot::geom_tile(color = "white")+
    ggplot::scale_fill_gradient2(low = "blue", high = "red",
                         midpoint = median_corr, limit = c(min_corr,1), space = "Lab",
                         name = paste(stringr::str_to_title(method), "\nCorrelation", sep = " ")) +
    ggplot::theme_minimal() +
    ggplot::theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
    ggplot::coord_fixed()
}

corr_scatter <- function(df, sample_1, sample_2, method = "spearman") {
  ggpubr::ggscatter(df, x = sample_1, y = sample_2,
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = method)
}
