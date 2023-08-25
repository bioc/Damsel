#' Drosophila melanogaster GATC region file
#'
#' Position information of the GATC regions across the dm6 genome. Regions are consecutive as the start and ends are part of the GATC sites.
#'
#' @format ## `regions_gatc_drosophila_dm6`
#' A data frame with 373,710 rows and 6 columns:
#' \describe{
#'   \item{Position}{Chromosome name and start of region}
#'   \item{seqnames}{Chromosome name}
#'   \item{start, end, width}{start and end position of region, width of the region}
#'   \item{strand}{Strand region is on}
#'   ...
#' }
#' @source <https://owenjm.github.io/damidseq_pipeline/>
"regions_gatc_drosophila_dm6"
