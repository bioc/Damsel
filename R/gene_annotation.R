
#' Extract gene information from Ensembl using biomaRt
#'
#' @description
#'`get_biomart_genes` accesses the Ensembl database via [biomaRt::useEnsembl() and biomaRt::getBM()] to obtain the location of genes from the selected species.
#' * Also identifies the number of GATC regions matching to each gene.
#'
#' @param species The species of interest. Format is first letter of genus, followed by full name of species, followed by gene_ensembl. For example: Drosophila melanogaster is dmelanogaster_gene_ensembl
#' @param version The Ensembl version of the genome. If not specified, default is 109 (the most recent update to the Drosophila melanogaster dm6 genome)
#' @param regions A data frame of GATC regions. If not specified, default is GATC regions from Drosophila melanogaster - dm6.
#'
#' @return A data.frame of information about the genes. Columns include: seqnames, start, end, width, strand, ensembl_gene_id, gene_name, ensembl_transcript_id, TSS (transcription start site), n_regions (number of overlapping GATC regions)
#' @export
#'
#' @examples
#' genes <- get_biomart_genes(species = "dmelanogaster_gene_ensembl",
#'                            version = 109,
#'                            regions = regions_gatc_drosophila_dm6)
#' head(genes)
#getGenes
get_biomart_genes <- function(species, version=109, regions=regions_gatc_drosophila_dm6) {
  if(!is.character(species)) {
    stop("Species must be a character vector")
  }
  if(!is.data.frame(regions)) {
    stop("Regions must be a data frame")
  }
  if(missing(version)) {
    message("Default version 109 used")
  }
  if(missing(regions)) {
    message("Default of drosophila dm6 regions used")
  }
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = species, version = version)
  BM.info <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical"),
                            filters = "chromosome_name",
                            values = unique(regions$seqnames),
                            mart = ensembl)
  gene_features <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
                                                 "chromosome_name", "start_position", "end_position", "strand",
                                                 "transcription_start_site"),
                                  filters = "ensembl_transcript_id",
                                  values = dplyr::filter(BM.info, .data$transcript_is_canonical == 1)$ensembl_transcript_id,
                                  mart = ensembl)
  gene_features <- gene_features %>% .[order(.$chromosome_name, .$start_position),]
  colnames(gene_features) <- c("ensembl_gene_id", "gene_name", "ensembl_transcript_id", "seqnames", "start", "end", "strand", "TSS")
  gene_features <- gene_features %>%
      plyranges::as_granges() %>%
      data.frame()

  overlap <- plyranges::find_overlaps_within(plyranges::as_granges(regions), plyranges::as_granges(gene_features)) %>%
      data.frame() %>%
      dplyr::group_by(.data$ensembl_gene_id) %>%
      dplyr::summarise(n_regions = dplyr::n()) %>%
      data.frame()
  gene_features$n_regions <- overlap[match(gene_features$ensembl_gene_id, overlap$ensembl_gene_id), "n_regions"]
  gene_features <- gene_features %>%
      dplyr::mutate(n_regions = dplyr::coalesce(.data$n_regions, 0),
                    seqnames = paste0("chr", .data$seqnames))

  gene_features
}


#' OLD: Annotating peaks to their closest gene
#'
#' `gene_annotate` takes both the peaks and genes as input and returns the paired results.
#' * The output of this function is a very large data frame with some confusing columns, the following function will help simplify the results.
#'
#' @param peaks data.frame of peaks as outputted from [aggregate_peaks]
#' @param genes data frame of gene information as outputted from [get_biomart_genes]
#'
#' @return data.frame of annotated peaks and genes.
#' Contains all columns of output of get_biomart_genes (with gene added as a prefix to seqnames, start, end, width, strand),
#' and all columns of output from aggregate_peaks.
#' Contains additional columns providing the: peak midpoint, distance from the midpoint to the start of the gene, etc
#' @export
gene_annotate <- function(peaks, genes) {
  if(!is.data.frame(peaks)) {
    stop("Require data.frame of peaks as outputted from `aggregate_peaks")
  }
  if(!is.data.frame(genes)) {
    stop("Requires data.frame of genes as outputted from `get_biomart_genes")
  }
  annotated <- plyranges::pair_nearest(plyranges::as_granges(genes), plyranges::as_granges(peaks)) %>%
      data.frame() %>%
      dplyr::mutate(peak_midpoint = (granges.y.start + granges.y.end)/2,
                    distance_to_start = granges.x.start - peak_midpoint,
                    words = ifelse(distance_to_start >= 0, "Upstream", "Downstream"),
                    combo = paste0(ensembl_gene_id, "-", number))
  overlaps_only <- plyranges::pair_overlaps(plyranges::as_granges(genes), plyranges::as_granges(peaks)) %>%
      data.frame() %>%
      dplyr::mutate(combo = paste0(ensembl_gene_id, "-", number))
  annotated <- annotated %>%
      dplyr::mutate(location = ifelse(combo %in% overlaps_only$combo, "Overlap", "NA"),
                    distance = dplyr::case_when(location == "NA" & words == "Upstream" ~ granges.x.start - granges.y.end,
                                         location == "NA" & words == "Downstream" ~ granges.y.start - granges.x.end,
                                         TRUE ~ 0))
  annotated <- annotated %>%
      dplyr::filter(distance <= 5000) %>%
      dplyr::group_by(number) %>%
      dplyr::mutate(peaks = dplyr::n(),
                    min_distance = ifelse(peaks == 1, distance_to_start, min(abs(distance_to_start))))
  annotated
}


#' OLD: Tabular display of peak statistical information and their closest genes
#'
#' `gene_annotate_organised` simplifies the output from `gene_annotate` to a clean format.
#'
#' @param annotated_peaks data.frame of annotated peaks [gene_annotate()]
#'
#' @return organised data frame of results with statistical information about each peak, alongside a list of the closest genes.
#' * contains information about the peak and gene location (for the closest gene), and provides a list of other nearby genes and their distance to the peak.
#' @export
gene_annotate_organised <- function(annotated_peaks) {
  if(!is.data.frame(annotated_peaks)) {
    stop("Requires data.frame of annotated peaks as outputted from `gene_annotate")
  }
  gene_peak <- annotated_peaks
  gene_peak <- gene_peak %>%
      dplyr::group_by(consec_dm) %>%
      dplyr::mutate(is_fbgn = stringr::str_detect(ensembl_gene_id, "FBgn"),
             n_genes = dplyr::n(),
             is_0 = ifelse(distance == 0, TRUE, FALSE),
             min_distance = ifelse(abs(min_distance) == abs(distance_to_start), abs(min_distance), 0)) %>%
      dplyr::mutate(has_min = ifelse(min_distance > 0, TRUE, FALSE)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(consec_dm, is_fbgn) %>%
      dplyr::mutate(n_coding = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(ratio = n_coding/n_genes) %>%
      dplyr::group_by(consec_dm, is_0) %>%
      dplyr::mutate(n_0 = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(is_fbgn == TRUE | ratio == 1 | has_min == 1) %>%
      dplyr::mutate(closest = ifelse(has_min == TRUE, gene_name, NA))
  closest <- gene_peak %>%
      dplyr::filter(!is.na(closest)) %>%
      dplyr::distinct(consec_dm, closest) %>%
      dplyr::group_by(consec_dm) %>%
      dplyr::mutate(number = 1:dplyr::n()) %>%
      dplyr::filter(number == 1) %>%
      dplyr::ungroup() %>%
      data.frame()
  others <- gene_peak %>%
      dplyr::group_by(consec_dm, closest) %>%
      dplyr::mutate(num = 1:dplyr::n(), check = ifelse(closest == "Yes" & num == 1, 1, 0)) %>%
      dplyr::filter(check != 1) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(consec_dm) %>%
      dplyr::summarise(list_gene = toString(gene_name),
                       list_gene_id = toString(ensembl_gene_id),
                       distance = toString(distance)) %>%
      dplyr::ungroup() %>%
      data.frame()

  closest$list_gene <-  others[match(closest$consec_dm, others$consec_dm), "list_gene"]
  closest$list_gene_id <- others[match(closest$consec_dm, others$consec_dm), "list_gene_id"]
  closest$distance <- others[match(closest$consec_dm, others$consec_dm), "distance"]
  closest
}

#' New: annotation of peaks and genes
#'
#' @description
#' `annotate_genes` identifies the closest gene(s) for the peaks outputted from `aggregate_peaks()`.
#'
#' This distance is relative, as the function will identify the closest genes, even if they are up to a million bp away. The max_distance parameter limits this, with a default setting of 5000 bp. All of the possible pairings are visible with `max_distance=NULL`.
#' The minimum distance between the peak and gene is calculated, (0 if the peak is within the gene or vice versa) and the relative position of the peak to the gene is also provided (Upstream, Downstream, Overlapping upstream, Contained within etc).
#'
#' @param peaks A data.frame of peaks as outputted from [aggregate_peaks()]
#' @param genes A data.frame of genes as outputted from [get_biomart_genes()]
#' @param regions A data.frame of GATC regions. If not specified, default is [regions_gatc_drosophila_dm6]
#' @param max_distance A number providing the limit for the minimum distance from peak to gene.
#' * Default is 5000. If set to `NULL`, will output all available combinations.
#'
#' @return A `list` of 3 `data.frames`:
#' * closest - every peak with it's closest gene
#' * top_5 - every peak with list of 5 closest genes
#' * all - all genes matching to each peak and all information
#' @export
#' @examples
#' set.seed(123)
#' example_regions <- random_regions()
#' dm_results <- random_edgeR_results()
#' peaks <- aggregate_peaks(dm_results)
#'
#' genes <- get_biomart_genes(species = "dmelanogaster_gene_ensembl",
#'                            version = 109,
#'                            regions = example_regions)
#'
#' annotate_genes(peaks, genes, example_regions, max_distance = 5000)
#' #view all combinations
#' annotate_genes(peaks, genes, example_regions, max_distance = NULL)
annotate_genes <- function(peaks, genes, regions=regions_gatc_drosophila_dm6, max_distance=5000) {
  if(!is.data.frame(peaks)) {
    stop("Require data.frame of peaks as outputted from `aggregate_peaks")
  }
  if(!is.data.frame(genes)) {
    stop("Requires data.frame of genes as outputted from `get_biomart_genes")
  }
  if(!is.data.frame(regions)) {
    stop("Requires data.frame of regions")
  }
  if(missing(regions)) {
    message("No regions provided, default regions_gatc_drosophila_dm6 used instead")
  }
  genes_gr <- plyranges::as_granges(genes)
  peaks_gr <- plyranges::as_granges(peaks)

  pair_gene <- plyranges::pair_nearest(genes_gr, peaks_gr) %>%
    data.frame()
  pair_peak <- plyranges::pair_nearest(peaks_gr, genes_gr) %>%
    data.frame()

  pair_gene <- pair_gene %>%
    stats::setNames(c("gene_seqnames", "gene_start", "gene_end", "gene_width", "gene_strand",
               "peak_seqnames", "peak_start", "peak_end", "peak_width", "peak_strand",
               colnames(.[11:ncol(.)])))
  pair_peak <- pair_peak %>%
    stats::setNames(c("peak_seqnames", "peak_start", "peak_end", "peak_width", "peak_strand",
               "gene_seqnames", "gene_start", "gene_end", "gene_width", "gene_strand",
               colnames(.[11:ncol(.)])))


  pair_peak <- pair_peak[,colnames(pair_gene)]

  pair_gene$from <- "gene"
  pair_peak$from <- "peak"

  overlaps <- plyranges::pair_overlaps(genes_gr, peaks_gr) %>%
    data.frame()
  overlaps <- overlaps %>%
    stats::setNames(c("gene_seqnames", "gene_start", "gene_end", "gene_width", "gene_strand",
               "peak_seqnames", "peak_start", "peak_end", "peak_width", "peak_strand",
               colnames(.[11:ncol(.)])))

  overlaps$from <- "overlap"
  overlaps$peak_order <- stringr::str_extract(overlaps$peak_id , "(?<=_)[0-9]+") %>%
    as.double()

  combo <- rbind(pair_gene, pair_peak)
  combo$peak_order <- stringr::str_extract(combo$peak_id , "(?<=_)[0-9]+") %>%
    as.double()
  combo <- combo %>%
    .[order(.$peak_order),]

  combo <- combo %>%
    dplyr::group_by(.data$peak_id, .data$ensembl_gene_id) %>%
    dplyr::mutate(count = dplyr::n(),
                  from = ifelse(.data$count == 1, .data$from, "both")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.)

  combo <- combo[,colnames(overlaps)]
  combo <- rbind(combo, overlaps) %>%
    data.frame() %>%
    .[order(.$peak_order),]

  combo <- combo %>%
    dplyr::group_by(.data$peak_id, .data$ensembl_gene_id) %>%
    dplyr::mutate(count = dplyr::n(),
                  from = ifelse(.data$count == 1, .data$from, "overlap")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.) %>%
    .[,1:24] %>%
    dplyr::mutate(seqnames = .data$gene_seqnames,
                  start = ifelse(.data$gene_start < .data$peak_start, .data$gene_start, .data$peak_start),
                  end = ifelse(.data$gene_end > .data$peak_end, .data$gene_end, .data$peak_end),
                  width = .data$end - .data$start + 1) %>%
    dplyr::mutate(peak_midpoint = (.data$peak_start + .data$peak_end)/2,
                  distance_TSS = .data$TSS - .data$peak_midpoint,
                  midpoint_is = ifelse(.data$distance_TSS >= 0, "Upstream", "Downstream")) %>%
    dplyr::mutate(abs_TSS = abs(.data$distance_TSS)) %>%
    dplyr::group_by(.data$peak_id) %>%
    dplyr::mutate(n_genes = dplyr::n()) %>%
    dplyr::ungroup()

  combo <- combo %>%
    dplyr::mutate(gap_ups = .data$gene_start - .data$peak_end,
                  gap_dow = .data$peak_start - .data$gene_end,
                  gap_st = .data$gene_start - .data$peak_start,
                  gap_en = .data$peak_end - .data$gene_end) %>%
    dplyr::mutate(position = dplyr::case_when(.data$gap_ups > 0 & .data$gap_ups <= .data$gap_st & .data$gap_dow < 0 & .data$gap_en < 0 ~ "Peak_upstream",
                                       .data$gap_ups < 0 & .data$gap_st < 0 & .data$gap_dow > 0 & .data$gap_dow <= .data$gap_en ~ "Peak_downstream",
                                       .data$gap_ups < 0 & .data$gap_st > 0 & .data$gap_dow < 0 & .data$gap_en < 0 ~ "Peak_overlap_upstream",
                                       .data$gap_ups < 0 & .data$gap_st > 0 & .data$gap_dow < 0 & .data$gap_en > 0 ~ "Peak_encompass_gene",
                                       .data$gap_ups < 0 & .data$gap_st < 0 & .data$gap_dow < 0 & .data$gap_en > 0 ~ "Peak_overlap_downstream",
                                       .data$gap_ups < 0 & .data$gap_st < 0 & .data$gap_dow < 0 & .data$gap_en < 0 ~ "Peak_within_gene"),
                  min_distance = dplyr::case_when(.data$position == "Peak_upstream" ~ .data$gap_ups,
                                                  .data$position == "Peak_downstream" ~ .data$gap_dow,
                                                  .data$position == "Peak_overlap_upstream" ~ .data$gap_st,
                                                  .data$position == "Peak_overlap_downstream" ~ .data$gap_en,
                                                  TRUE ~ 0)) %>%
    .[,!(colnames(.) %in% c("gap_ups", "gap_dow", "gap_st", "gap_en"))]

  combo <- combo %>%
    dplyr::mutate(num = 1:dplyr::n())

  n_regions <- combo[,c("seqnames", "start", "end", "width")] %>%
    dplyr::mutate(num = 1:dplyr::n())

  n_regions <- plyranges::find_overlaps_within(plyranges::as_granges(dplyr::mutate(regions,
                                                                                   seqnames = paste0("chr", .data$seqnames))),
                                               plyranges::as_granges(n_regions)) %>%
    data.frame() %>%
    dplyr::group_by(.data$num) %>%
    dplyr::summarise(n_region = dplyr::n()) %>%
    data.frame()
  combo$total_regions <- n_regions[match(combo$num, n_regions$num), "n_region"]

  combo <- combo %>%
    .[order(.$peak_order, .$min_distance),] %>%
    dplyr::group_by(.data$peak_id) %>%
    dplyr::mutate(count = 1:dplyr::n()) %>%
    dplyr::ungroup()

  if(!is.null(max_distance)) {
    combo <- dplyr::filter(combo, .data$min_distance <= max_distance)
  }

  col_order <- c("seqnames", "start", "end", "width",
                 "total_regions", "n_regions_dm", "peak_id", "rank_p",
                 "gene_position", "ensembl_gene_id", "gene_name", "midpoint_is",
                 "position")

  closest <- combo %>%
    dplyr::filter(.data$count == 1) %>%
    dplyr::mutate(gene_position = paste0(.data$gene_seqnames, ":", .data$gene_strand, ":", .data$gene_start, "-", .data$gene_end),
                  peak_position = paste0(.data$peak_seqnames, ":", .data$peak_start, "-", .data$peak_end)) %>%
    .[,col_order]

  combo_list <- combo %>%
    dplyr::filter(.data$count <= 5) %>%
    dplyr::group_by(.data$peak_id) %>%
    dplyr::mutate(start_region = min(.data$gene_start),
                  end_region = max(.data$gene_end),
                  start = ifelse(.data$start_region < .data$start, .data$start_region, .data$start),
                  end = ifelse(.data$end_region > .data$end, .data$end_region, .data$end)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$seqnames, .data$start, .data$end, .data$peak_id, .data$rank_p, .data$n_genes) %>%
    dplyr::summarise(list_ensembl = toString(.data$ensembl_gene_id),
                     list_gene = toString(.data$gene_name),
                     position = toString(.data$position),
                     distance_TSS = toString(.data$distance_TSS),
                     min_distance = toString(.data$min_distance)) %>%
    dplyr::ungroup()
    cols_start <- which(colnames(combo) %in% c("seqnames", "start", "end", "width"))

    combo <- combo[,c(cols_start, ncol(combo),
                      1:(utils::head(cols_start, 1)-1),
                      (utils::tail(cols_start, 1)+1):(ncol(combo)-1))] %>%
      .[,!(colnames(.) %in% c("from", "peak_order", "abs_TSS", "count", "num"))]

  list_results <- list(closest = closest, top_five = combo_list, all = combo)
  list_results
}




extend_upstream <- function(granges, upstream=0) {
  on_plus <- strand(granges) == "+"
  new_start <- start(granges) - ifelse(on_plus, upstream, 0)
  new_end <- end(granges) + ifelse(on_plus, 0, upstream)
  start(granges) <- new_start
  end(granges) <- new_end
  granges
}


