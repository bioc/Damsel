#' Obtain region counts for BAM files
#'
#' `countBamInGATC()` obtains the raw counts for the regions between GATC sites, from indexed BAM files specified in the path.
#'
#' @param path_to_bams A string identifying the directory containing the BAM files.
#' @param regions A GRanges object of GATC regions. The GATC regions can be made with `gatc_track()`.
#' @param nthreads The number of computer cores to be used to parallelise the function and decrease it's run time. If not specified, will use default (2 cores).
#' * If computer is being used for multiple tasks at once, we recommend reducing the number of cores - or leave it at the default setting.
#' * The number of available cores can be checked using [parallel::detectCores()]
#' @param ... Other arguments passed onto `Rsubread::featureCounts()`
#'
#' @return A `data.frame` containing the GATC region information in the form in the columns: seqnames (chromosome), start, end, width, and strand. The count information for the BAM files is in the subsequent columns, named by the name of the BAM file.
#' * The ".bam" extension is retained in the sample name as an identifier for the sample columns
#' * If necessary, at this stage please rearrange the BAM file columns so they are ordered in the following way: Dam_1, Fusion_1, Dam_2, Fusion_2 etc
#' * The DamID data captures the ~75bp region extending from each GATC site, so although regions are of differing widths, there is a null to minimal length bias present on the data, and does not require length correction.
#' @references Liao Y, Smyth GK, Shi W (2019). “The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads.” Nucleic Acids Research, 47, e47. doi:10.1093/nar/gkz114.
#' Morgan M, Pagès H, Obenchain V, Hayden N (2024). Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import. R package version 2.19.3, https://bioconductor.org/packages/Rsamtools.
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' example_regions <- random_regions()
#' counts.df <- countBamInGATC(path_to_bams,
#'     regions = example_regions,
#'     nthreads = 2
#' )
#' head(counts.df)
#' # rearrange columns of bam files so that: Dam_1, Fusion_1, Dam_2, Fusion_2
#' @export
countBamInGATC <- function(path_to_bams, regions, nthreads=2, ...) {
    if (!is.character(path_to_bams)) {
        stop("Path to bams must be a character vector")
    }
    if (!is.data.frame(regions) && !inherits(regions, "GRanges")) {
        stop("GATC region file must be a GRanges object")
    }
    if (Sys.info()["sysname"] == "Windows") {
      nthreads <- 1
    }
    regions <- data.frame(regions)
    if(!"NCBI" %in% GenomeInfoDb::seqlevelsStyle(unfactor(unique(regions$seqnames)))) {
        regions <- ..changeStyle(regions, "NCBI")
    }
    if(!"Position" %in% colnames(regions)) {
        regions <- regions %>% dplyr::mutate(Position = paste0("chr", seqnames, "-", start))
    }
    regions <- regions[,c("Position", "seqnames", "start", "end", "width", "strand"), drop = FALSE]
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
    path_to_bams <- ifelse(grepl("/$", path_to_bams),
        gsub("/$", "", path_to_bams), path_to_bams
    )
    list_files <- file.path(path_to_bams, list_files)

    scan_result <- ..checkPaired(list_files)
    paired_samples <- list_files[grep(TRUE, scan_result)]
    single_samples <- list_files[grep(FALSE, scan_result)]

    counts_feature <- regions
    names_chrom <- names(Rsamtools::scanBamHeader(list_files[1])[[1]][[1]])
    same_name <- counts_feature$seqnames %in% names_chrom %>% unique()

    if (same_name == FALSE) {
        counts_feature$seqnames <- paste0("chr", counts_feature$seqnames)
    }
    if (length(paired_samples) != 0) {
        counts_feature <- cbind(
            counts_feature,
            Rsubread::featureCounts(paired_samples,
                annot.ext = regions_feat,
                isPairedEnd = TRUE, allowMultiOverlap = TRUE, fraction = TRUE,
                nthreads = nthreads, ...
            )$counts
        )
    }
    if (length(single_samples) != 0) {
        counts_feature <- cbind(
            counts_feature,
            Rsubread::featureCounts(single_samples,
                annot.ext = regions_feat,
                allowMultiOverlap = TRUE, fraction = TRUE, nthreads = nthreads,
                ...
            )$counts
        )
    }
    if (same_name == TRUE) {
        counts_feature$seqnames <- paste0("chr", counts_feature$seqnames)
    }

    counts_feature
}

..checkPaired <- function(list_files) {
    scan_result <- vapply(list_files, Rsamtools::testPairedEndBam, c(a = TRUE))
    scan_result
}
