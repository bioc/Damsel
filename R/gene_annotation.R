#getGenes
get_biomart_genes <- function(species, version=109, regions) {
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = species, version = version)
  BM.info <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical"),
                            filters = "chromosome_name",
                            values = unique(regions$seqnames),
                            mart = ensembl)
  gene_features <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
                                                 "chromosome_name", "start_position", "end_position", "strand",
                                                 "transcription_start_site"),
                                  filters = "ensembl_transcript_id",
                                  values = filter(BM.info, transcript_is_canonical == 1)$ensembl_transcript_id,
                                  mart = ensembl)
  gene_features <- gene_features %>% .[order(.$chromosome_name, .$start_position),]
  colnames(gene_features) <- c("ensembl_gene_id", "gene_name", "ensembl_transcript_id", "seqnames", "start", "end", "strand", "TSS")
  gene_features <- gene_features %>%
                     plyranges::as_granges() %>%
                     data.frame()

  overlap <- plyranges::find_overlaps_within(plyranges::as_granges(regions), plyranges::as_granges(gene_features)) %>%
               data.frame() %>%
               dplyr::group_by(ensembl_gene_id) %>%
               dplyr::summarise(n_regions = n()) %>%
               data.frame()
  gene_features$n_regions <- overlap[match(gene_features$ensembl_gene_id, overlap$ensembl_gene_id), "n_regions"]
  gene_features <- gene_features %>% dplyr::mutate(n_regions = dplyr::coalesce(n_regions, 0), seqnames = paste0("chr", seqnames))

  gene_features
}



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
