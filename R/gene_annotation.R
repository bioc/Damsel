
#' Get gene information from Ensembl using biomaRt
#'
#'`get_biomart_genes` accesses the Ensembl database via [biomaRt::useEnsembl() and biomaRt::getBM()] to obtain the location of genes from the selected species.
#' Also identifies the number of GATC regions matching to each gene.
#'
#' @param species species of interest. Format is first letter of genus, followed by full name of species, followed by gene_ensembl. For example: Drosophila melanogaster is dmelanogaster_gene_ensembl
#' @param version Ensembl version of genome. Default is 109 (the most recent update to the Drosophila melanogaster dm6 genome)
#' @param regions data frame of GATC regions. Default is GATC regions from Drosophila melanogaster - dm6.
#'
#' @return data.frame of information about the genes. Columns include (I don't know what order right now): chromosome, start, end, width, strand, gene name, ensembl_gene_id, ensembl_transcript_id, n_overlapping_gatc_regions
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
                                  values = dplyr::filter(BM.info, transcript_is_canonical == 1)$ensembl_transcript_id,
                                  mart = ensembl)
  gene_features <- gene_features %>% .[order(.$chromosome_name, .$start_position),]
  colnames(gene_features) <- c("ensembl_gene_id", "gene_name", "ensembl_transcript_id", "seqnames", "start", "end", "strand", "TSS")
  gene_features <- gene_features %>%
      plyranges::as_granges() %>%
      data.frame()

  overlap <- plyranges::find_overlaps_within(plyranges::as_granges(regions), plyranges::as_granges(gene_features)) %>%
      data.frame() %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(n_regions = dplyr::n()) %>%
      data.frame()
  gene_features$n_regions <- overlap[match(gene_features$ensembl_gene_id, overlap$ensembl_gene_id), "n_regions"]
  gene_features <- gene_features %>%
      dplyr::mutate(n_regions = dplyr::coalesce(n_regions, 0),
                    seqnames = paste0("chr", seqnames))

  gene_features
}


#' Annotating peaks to their closest gene
#'
#' `gene_annotate` takes both the peaks and genes as input and returns the paired results.
#'
#' @param peaks data.frame of peaks as outputted from [aggregate_peaks]
#' @param genes data frame of gene information as outputted from [get_biomart_genes]
#'
#' @return data.frame of annotated peaks and genes in a horrible messy format I don't like to look at.
#' @export
#'
#' @examples
#' # set up peaks
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' dge <- edgeR_set_up(counts.df)
#' de_results <- edgeR_results(dge, p.value = 0.05, lfc = 1)
#' peaks <- aggregate_peaks(de_results, regions = regions_gatc_drosophila_dm6)
#' # set up genes
#' genes <- get_biomart_genes(species = "dmelanogaster_gene_ensembl",
#'                            version = 109,
#'                            regions = regions_gatc_drosophila_dm6)
#'
#' annotated_peaks <- gene_annotate(peaks, genes)
#' head(annotated_peaks)
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


#' Tabular display of peak statistical information and their closest genes
#'
#' `gene_annotate_organised` simplifies the output from `gene_annotate` to a clean format. meaning I probably want it all as 1 fn
#'
#' @param annotated_peaks data.frame of annotated messy peaks I don't like from [gene_annotate()]
#'
#' @return organised data frame of results with statistical information about each peak, alongside a list of the closest genes.
#' @export
#'
#' @examples
#' # set up peaks
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' dge <- edgeR_set_up(counts.df)
#' de_results <- edgeR_results(dge, p.value = 0.05, lfc = 1)
#' peaks <- aggregate_peaks(de_results, regions = regions_gatc_drosophila_dm6)
#' # set up genes
#' genes <- get_biomart_genes(species = "dmelanogaster_gene_ensembl", version = 109, regions = regions_gatc_drosophila_dm6)
#' annotated_peaks <- gene_annotate(peaks, genes)
#'
#' annotated_peaks <- gene_annotate_organised(annotated_peaks)
#' head(annotated_peaks)
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
