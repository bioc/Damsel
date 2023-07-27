extend_upstream <- function(granges, upstream=0) {
  on_plus <- strand(granges) == "+"
  new_start <- start(granges) - ifelse(on_plus, upstream, 0)
  new_end <- end(granges) + ifelse(on_plus, 0, upstream)
  start(granges) <- new_start
  end(granges) <- new_end
  granges
}

gene_annotate_fn <- function(granges, upstream, peaks) {
  on_plus <- strand(granges) == "+"
  new_start <- start(granges) - ifelse(on_plus, upstream, 0)
  new_end <- end(granges) + ifelse(on_plus, 0, upstream)
  start(granges) <- new_start
  end(granges) <- new_end

  plyranges::find_overlaps_within(plyranges::as_granges(dplyr::mutate(peaks, start = peak_start, end = peak_end)),
                       plyranges::as_granges(dplyr::mutate(data.frame(granges), gene_width = width))) %>%
    data.frame()

}
