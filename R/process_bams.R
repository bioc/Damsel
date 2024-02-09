#' Obtain region counts for BAM files
#'
#' `process_bams()` obtains the raw counts for the regions between GATC sites, from indexed BAM files specified in the path.
#'
#' @param path_to_bams A character string input vector. This must identify the directory containing the BAM files.
#' @param regions A data.frame of GATC regions. The GATC regions can be made with `gatc_track()`.
#' @param nthreads The number of computer cores to be used to parallelise the function and decrease it's run time. If not specified, will use default (2 cores).
#' * If computer is being used for multiple tasks at once, we recommend reducing the number of cores - or leave it at the default setting.
#' * The number of available cores can be checked using [parallel::detectCores()]
#' @param ... Other arguments passed onto `Rsubread::featureCounts()`
#'
#' @return A `data.frame` containing the GATC region information in the form in the columns: seqnames (chromosome), start, end, width, and strand. The count information for the BAM files is in the subsequent columns, named by the name of the BAM file.
#' * The ".bam" extension is retained in the sample name as an identifier for the sample columns
#' * If necessary, at this stage please rearrange the BAM file columns so they are ordered in the following way: Dam_1, Fusion_1, Dam_2, Fusion_2 etc
#' * The DamID data captures the ~75bp region extending from each GATC site, so although regions are of differing widths, there is a null to minimal length bias present on the data, and does not require length correction.
#'
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' example_regions <- random_regions()
#' counts.df <- process_bams(path_to_bams,
#'     regions = example_regions,
#'     nthreads = 2
#' )
#' head(counts.df)
#' # rearrange columns of bam files so that: Dam_1, Fusion_1, Dam_2, Fusion_2
#' head(counts.df[, c(1:6, 9, 7, 10, 8)])
#' @export
# rename to processBams
process_bams <- function(path_to_bams, regions, nthreads = 2, ...) {
    if (!requireNamespace("Rsubread", quietly = TRUE)) {
        stop("Package \"Rubsread\" must be installed to use this function.",
            call. = FALSE
        )
    }
    if (!is.character(path_to_bams)) {
        stop("Path to bams must be a character vector")
    }
    if (!is.data.frame(regions)) {
        stop("GATC region file must be a data.frame")
    }
    if (missing(nthreads)) {
        message("nthreads missing, 2 cores used instead")
    }

    regions_feat <- regions[, !names(regions) %in% "width"]
    colnames(regions_feat) <- c("GeneID", "Chr", "Start", "End", "Strand")
    regions_feat$Start <- as.integer(regions_feat$Start)
    regions_feat$End <- as.integer(regions_feat$End)

    list_files <- list.files(path_to_bams, pattern = ".bam", ignore.case = TRUE)
    # check
    if (length(list_files) == 0) {
        stop("No bam files present in path")
    }
    if (length(list_files[grepl(".bai", list_files)]) == 0) {
        stop("No .bai files present in path: Please ensure every BAM file has a corresponding .bai file")
    }
    list_files <- list_files[!grepl(".bai", list_files)]
    path_to_bams <- ifelse(substring(path_to_bams, first = nchar(path_to_bams)) == "/", path_to_bams, paste0(path_to_bams, "/"))
    list_files <- paste0(path_to_bams, list_files)

    scan_result <- check_paired(list_files)
    paired_samples <- list_files[grep(TRUE, scan_result)]
    single_samples <- list_files[grep(FALSE, scan_result)]

    counts_feature <- regions
    names_chrom <- names(Rsamtools::scanBamHeader(list_files[1])[[1]][[1]])
    same_name <- counts_feature$seqnames %in% names_chrom %>% unique()

    if (same_name == FALSE) {
        counts_feature$seqnames <- paste0("chr", counts_feature$seqnames)
    }
    if (length(paired_samples) != 0) {
        for (i in seq_len(length(paired_samples))) {
            counts_feature <- cbind(counts_feature, data.frame(Rsubread::featureCounts(paired_samples[i], annot.ext = regions_feat, isPairedEnd = TRUE, allowMultiOverlap = TRUE, fraction = TRUE, nthreads = nthreads, ...)$counts))
        }
    }
    if (length(single_samples) != 0) {
        for (i in seq_len(length(single_samples))) {
            counts_feature <- cbind(counts_feature, data.frame(Rsubread::featureCounts(single_samples[i], annot.ext = regions_feat, allowMultiOverlap = TRUE, fraction = TRUE, nthreads = nthreads, ...)$counts))
        }
    }
    if (same_name == TRUE) {
        counts_feature$seqnames <- paste0("chr", counts_feature$seqnames)
    }

    counts_feature
}

check_paired <- function(list_files) {
    scan_result <- list()
    for (i in seq_len(length(list_files))) {
        scan_result <- list(scan_result, Rsamtools::testPairedEndBam(list_files[i]))
    }
    scan_result <- unlist(scan_result)
    scan_result
}


