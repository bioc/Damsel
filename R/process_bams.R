#' Obtain region counts for BAM files
#'
#' `process_bams()` obtains the raw counts for the regions between GATC sites, from indexed BAM files specified in the path.
#'
#' @param path_to_bams A character string input vector. This must identify the directory containing the BAM files.
#' @param regions A data.frame of GATC regions. If not specified, will use default data file `regions_gatc_drosophila_dm6`. The GATC regions can be made with `gatc_track()`.
#' @param cores The number of computer cores to be used to parallelise the function and decrease it's run time. If not specified, will use default (2 cores).
#' * If computer is being used for multiple tasks at once, we recommend reducing the number of cores - or leave it at the default setting.
#' * The number of available cores can be checked using [parallel::detectCores()]
#'
#' @return A data.frame containing the GATC region information in the form in the columns: seqnames (chromosome), start, end, width, and strand. The count information for the BAM files is in the subsequent columns, named by the name of the BAM file.
#' * The ".bam" extension is retained in the sample name as an identifier for the sample columns
#' * If necessary, at this stage please rearrange the BAM file columns so they are ordered in the following way: Dam_1, Fusion_1, Dam_2, Fusion_2 etc
#' * The DamID data captures the ~75bp region extending from each GATC site, so although regions are of differing widths, there is a null to minimal length bias present on the data, and does not require length correction.
#'
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores=2)
#' head(counts.df)
#' # rearrange columns of bam files so that: Dam_1, Fusion_1, Dam_2, Fusion_2
#' head(counts.df[,c(1:6,7,10,8,11,9,12)])
#' @export
#rename to processBams
process_bams <- function(path_to_bams, regions=regions_gatc_drosophila_dm6, cores=2) {
  if(!is.character(path_to_bams)) {
    stop("Path to bams must be a character vector")
  }
  if(!is.data.frame(regions)) {
    stop("GATC region file must be a data.frame")
  }
  if(missing(regions)) {
    message("regions missing, regions_gatc_drosophila_dm6 used instead")
  }
  if(missing(cores)) {
    message("cores missing, 2 cores used instead")
  }

  #list of all bam files
  files <- list.files(path_to_bams, pattern = ".bam")
  #check
  if(length(files) == 0) {
    stop("No bam files present in path")
  }

  files <- files %>%
      data.frame() %>%
      stats::setNames(., "file") %>%
      dplyr::mutate(bai = stringr::str_detect(.data$file, "bai"))
  if(S4Vectors::isEmpty(dplyr::filter(files, .data$bai == TRUE))) {
    stop("No .bai files present in path")
  }
  #want to add a check for unequal n of bams and bai

    #remove any bai files
  files <- files %>%
      dplyr::filter(.data$bai == FALSE)
  #add / to get path for each file
  path_to_bams <- ifelse(substring(path_to_bams, first = nchar(path_to_bams)) == "/", path_to_bams, paste0(path_to_bams, "/"))
  files <- paste0(path_to_bams, files$file)
  df <- regions
  #test format of seqnames and separate them
  names_chrom <- names(Rsamtools::scanBamHeader(files[1])[[1]][[1]])
  same_name <- df$seqnames %in% names_chrom %>% unique()
  if(same_name == TRUE) {
    one_count <- parallel::mclapply(files, function(x) {exomeCopy::countBamInGRanges(x, plyranges::as_granges(df))}, mc.cores = cores) %>%
      data.frame()
    colnames(one_count) <- files
  } else {
    df$seqnames <- paste0("chr", df$seqnames)
    one_count <- parallel::mclapply(files, function(x) {exomeCopy::countBamInGRanges(x, plyranges::as_granges(df))}, mc.cores = cores) %>%
      data.frame()
    colnames(one_count) <- files
  }

  df2 <- cbind(df, one_count)
  df2$seqnames <- paste0("chr", df2$seqnames)
  colnames(df2) <- sub(path_to_bams, "", colnames(df2))
  df2

}



