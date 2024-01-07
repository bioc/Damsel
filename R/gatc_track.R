
#' Extract GATC regions
#'
#' `gatc_region_fn` identifies and extracts the GATC sites and regions from a BSgenome object or a fasta file.
#'
#' @param object A BSgenome package OR the path to a FASTA file.
#'
#' @return A `list` object composed of two `data.frames`.
#'
#' The first; "regions" provides the GATC regions used in downstream analysis. The second, "sites", provides the positions of the GATC sites, and is used in plotting the results.
#' @export
#' @examples
#' if (require("BSgenome.Dmelanogaster.UCSC.dm6")) {
#'   gatc <- gatc_region_fn(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)
#'
#'   head(gatc$regions)
#'
#'   head(gatc$sites)
#' }
gatc_region_fn <- function(object) {
  if(missing(object)) {
    stop("Must have a BSgenome object such as BSgenome.Dmelanogaster.UCSC.dm6, OR the path to a FASTA file")
  }
  if(!class(object) %in% c("BSgenome", "character")) {
    stop("Must have a BSgenome object such as BSgenome.Dmelanogaster.UCSC.dm6, OR the path to a FASTA file")
  }
  if(class(object) %in% "BSgenome") {
    fasta <- object
    names_fasta <- GenomeInfoDb::seqnames(object)
  } else if(class(object) %in% "character") {
    fasta <- Biostrings::readDNAStringSet(object)
    names_fasta <- names(fasta)
  }

  length_names <- stringr::str_detect(names_fasta, "Y") %>%
    grep("TRUE", .,) %>%
    .[1]
  seq_names <- names_fasta[seq_len(length_names)]

  df <- lapply(seq_names, function(x) cbind(data.frame(Biostrings::matchPattern("GATC", fasta[[x]])),
                                 seqnames = x)) %>%
    dplyr::bind_rows()

  df <- df[,c("seqnames", "start", "end", "width")]
  df$seqnames <- gsub("chr", "", as.character(df$seqnames))
  df$strand <- "*"

  regions <- df %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::mutate(start = .data$start + 1,
                  end = dplyr::lead(.data$start)) %>%
    dplyr::filter(!is.na(.data$end)) %>%
    dplyr::mutate(width = .data$end - .data$start + 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Position = paste0("chr", seqnames, "-", start))
  regions <- regions[,c("Position", "seqnames", "start", "end", "width", "strand")]

  list(regions = regions, sites = df)

}
