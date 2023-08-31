#' Initial processing of BAM files
#'
#' `process_bams()` obtains the counts for the regions between GATC sites, from BAM files specified in the path.
#'
#' @param path_to_bams Character string input vector. Should identify directory containing the BAM files.
#' @param regions Data frame of GATC regions. Default is GATC regions from Drosophila melanogaster - dm6.
#' @param cores Specify computer cores to be used to parallelise function and speed output. Default is to use all available cores. If computer is being used for multiple tasks at once, reccommend reducing the number of cores. Number of available cores can be checked using [parallel::detectCores()]
#'
#' @return A data frame containing the GATC region information in the form in the columns: seqnames (chromosome), start, end, width, and strand. The count information for the BAM files is in the subsequent columns, named by the name of the BAM file.
#' * If necessary, at this stage please rearrange the BAM file columns so they are ordered in the following way: Dam_1, Fusion_1, Dam_2, Fusion_2 etc
#' * The DamID data captures the ~75bp region extending from each GATC site, so although regions are of differing widths, there is a null to minimal length bias present on the data, and does not require length correction.
#'
#' @examples
#'
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
    message("cores missing, [parallel::detectCores()] used instead")
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
      dplyr::mutate(bai = stringr::str_detect(.$file, "bai"))
  if(S4Vectors::isEmpty(dplyr::filter(files, bai == TRUE))) {
    stop("No .bai files present in path")
  }
  #want to add a check for unequal n of bams and bai

    #remove any bai files
  files <- files %>%
      dplyr::filter(bai == FALSE)
  #add / to get path for each file
  path_to_bams <- ifelse(substring(path_to_bams, first = nchar(path_to_bams)) == "/", path_to_bams, paste0(path_to_bams, "/"))
  files <- paste0(path_to_bams, files$file)
  df <- regions
  #test format of seqnames and separate them
  one_count <- parallel::mclapply(files, function(x) {exomeCopy::countBamInGRanges(x, plyranges::as_granges(df))}, mc.cores = cores) %>%
      data.frame()
  colnames(one_count) <- files
#want i in seq_len or something
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



