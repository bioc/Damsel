#' Deprecated: Extract gene information from Ensembl using biomaRt
#'
#' @description
#' This function has been deprecated and replaced with `collateGenes`
#' `get_biomart_genes` accesses the Ensembl database via [biomaRt::useEnsembl() and biomaRt::getBM()] to obtain the location of genes from the selected species.
#' * Also identifies the number of GATC regions matching to each gene.
#'
#' @param species The species of interest. Format is first letter of genus, followed by full name of species, followed by gene_ensembl. For example: Drosophila melanogaster is dmelanogaster_gene_ensembl
#' @param version The Ensembl version of the genome. If not specified, default is 109 (the most recent update to the Drosophila melanogaster dm6 genome)
#' @param regions A data frame of GATC regions.
#'
#' @return A data.frame of information about the genes. Columns include: seqnames, start, end, width, strand, ensembl_gene_id, gene_name, ensembl_transcript_id, TSS (transcription start site), n_regions (number of overlapping GATC regions)
#' @export
get_biomart_genes <- function(species, version = 109, regions) {
    .Deprecated("collateGenes")
    warning("This function has been deprecated in v0.7.0")
    if (!is.character(species)) {
        stop("Species must be a character vector")
    }
    if (!is.data.frame(regions)) {
        stop("Regions must be a data frame")
    }
    if (missing(version)) {
        message("Default version 109 used")
    }

    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = species, version = version)
    BM.info <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical"),
        filters = "chromosome_name",
        values = unique(regions$seqnames),
        mart = ensembl
    )
    gene_features <- biomaRt::getBM(
        attributes = c(
            "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
            "chromosome_name", "start_position", "end_position", "strand",
            "transcription_start_site"
        ),
        filters = "ensembl_transcript_id",
        values = dplyr::filter(BM.info, .data$transcript_is_canonical == 1)$ensembl_transcript_id,
        mart = ensembl
    )
    gene_features <- gene_features %>% .[order(.$chromosome_name, .$start_position), ]
    colnames(gene_features) <- c("ensembl_gene_id", "gene_name", "ensembl_transcript_id", "seqnames", "start", "end", "strand", "TSS")

    overlap <- plyranges::find_overlaps_within(plyranges::as_granges(regions), plyranges::as_granges(gene_features)) %>%
        data.frame() %>%
        dplyr::group_by(.data$ensembl_gene_id) %>%
        dplyr::summarise(n_regions = dplyr::n()) %>%
        data.frame()
    gene_features$n_regions <- overlap[match(gene_features$ensembl_gene_id, overlap$ensembl_gene_id), "n_regions"]
    gene_features <- gene_features %>%
        dplyr::mutate(
            n_regions = dplyr::coalesce(.data$n_regions, 0),
            seqnames = paste0("chr", .data$seqnames)
        )

    gene_features
}

#' New: Get list of genes
#'
#' Takes a Txdb object, path to a gff file, or a species (biomaRt) and returns a GRanges of genes.
#'
#' @param genes A Txdb object, path to file, or a species for accessing biomaRt.
#' @param regions GATC region file.
#' @param org.Db Required if using a Txdb object so to access gene names.
#' @param version Required for using biomaRt.
#'
#' @return A GRanges object of genes and available supplementary information - specifically the TSS, and number of GATC regions overlapping the gene.
#' @export
#' @references Carlson M (2019). org.Dm.eg.db: Genome wide annotation for Fly. R package version 3.8.2.
#' Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, 1184–1191.
#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A, Huber W (2005). “BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis.” Bioinformatics, 21, 3439–3440.
#' Lee, Stuart, Cook, Dianne, Lawrence, Michael (2019). “plyranges: a grammar of genomic data transformation.” Genome Biol., 20(1), 4. http://dx.doi.org/10.1186/s13059-018-1597-8.
#' Team BC, Maintainer BP (2019). TxDb.Dmelanogaster.UCSC.dm6.ensGene: Annotation package for TxDb object(s). R package version 3.4.6.
#' @examples
#' library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#' library(org.Dm.eg.db)
#' set.seed(123)
#' example_regions <- random_regions()
#' txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#' genes <- collateGenes(genes = txdb, regions = example_regions, org.Db = org.Dm.eg.db)
#' head(genes)
collateGenes <- function(genes, regions, org.Db=NULL, version=NULL) {
    regions <- data.frame(regions)
    if (inherits(genes, "TxDb")) {
        genes_ <- GenomicFeatures::genes(genes)
        genes_$ensembl_gene_id <- names(genes_)
        genes_ <- data.frame(genes_)[, c("seqnames", "start", "end", "width", "strand", "ensembl_gene_id")]
        gene_names <- AnnotationDbi::select(org.Db, keys = AnnotationDbi::keys(genes), keytype = "ENSEMBL", columns = c("SYMBOL"))
        genes_$gene_name <- gene_names[match(genes_$ensembl_gene_id, gene_names$ENSEMBL), "SYMBOL"]
        genes_ <- genes_ %>% dplyr::mutate(TSS = ifelse(.data$strand == "+", .data$start, .data$end))
        message("TSS taken as start of gene, taking strand into account")
        genes_ <- plyranges::as_granges(genes_)
        seq_chr <- GenomeInfoDb::seqlevels(genes_)[1] %>% grepl("chr", .)
        if (seq_chr == TRUE) {
            regions$seqnames <- paste0("chr", regions$seqnames)
        }
    } else if (grep("ensembl", genes)) {
        genes_ <- ..accessBiomart(species = genes, regions, version)
    }
    genes_ <- ..annotateGeneGatc(genes_, regions)
    genes_
}

..accessBiomart <- function(species, regions, version) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("Package \"biomaRt\" must be installed to use this function.",
            call. = FALSE
        )
    }
    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = species, version = version)
    BM.info <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical"),
        filters = "chromosome_name",
        values = unique(regions$seqnames),
        mart = ensembl
    )
    gene_features <- biomaRt::getBM(
        attributes = c(
            "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
            "chromosome_name", "start_position", "end_position", "strand",
            "transcription_start_site"
        ),
        filters = "ensembl_transcript_id",
        values = dplyr::filter(BM.info, .data$transcript_is_canonical == 1)$ensembl_transcript_id,
        mart = ensembl
    )
    gene_features <- gene_features %>% .[order(.$chromosome_name, .$start_position), ]
    colnames(gene_features) <- c("ensembl_gene_id", "gene_name", "ensembl_transcript_id", "seqnames", "start", "end", "strand", "TSS")
    gene_features <- plyranges::as_granges(gene_features)
}

..annotateGeneGatc <- function(genes, regions) {
    # regions$seqnames <- paste0("chr", regions$seqnames)
    overlap <- plyranges::find_overlaps_within(plyranges::as_granges(regions), plyranges::as_granges(genes)) %>%
        data.frame() %>%
        dplyr::group_by(.data$ensembl_gene_id) %>%
        dplyr::summarise(n_regions = dplyr::n()) %>%
        data.frame()
    genes_ <- data.frame(genes)
    genes_$n_regions <- overlap[match(genes_$ensembl_gene_id, overlap$ensembl_gene_id), "n_regions"]
    genes_ <- genes_ %>%
        dplyr::mutate(
            n_regions = dplyr::coalesce(.data$n_regions, 0)
            # seqnames = ifelse(grepl("chr", seqnames), seqnames, paste0("chr", seqnames))
        ) %>%
        .[order(.$start), ] %>%
        .[order(.$seqnames), ]
    genes_ <- plyranges::as_granges(genes_)
    genes_
}

#' Annotation of peaks and genes
#'
#' @description
#' `annotate_genes` identifies the closest gene(s) for the peaks outputted from `aggregate_peaks()`.
#'
#' This distance is relative, as the function will identify the closest genes, even if they are up to a million bp away. The max_distance parameter limits this, with a default setting of 5000 bp. All of the possible pairings are visible with `max_distance=NULL`.
#' The minimum distance between the peak and gene is calculated, (0 if the peak is within the gene or vice versa) and the relative position of the peak to the gene is also provided (Upstream, Downstream, Overlapping upstream, Contained within etc).
#'
#' @param peaks A data.frame of peaks as outputted from [aggregate_peaks()].
#' @param genes A data.frame of genes as outputted from [get_biomart_genes()].
#' @param regions A `GRanges` object of GATC regions.
#' @param max_distance A number providing the limit for the minimum distance from peak to gene.
#' * Default is 5000. If set to `NULL`, will output all available combinations.
#'
#' @return A `list` of 3 `data.frames`:
#' * closest - every peak with it's closest gene
#' * top_5 - every peak with list of 5 closest genes
#' * all - all genes matching to each peak and all information
#' @export
#' @examples
#' library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#' library(org.Dm.eg.db)
#' set.seed(123)
#' example_regions <- random_regions()
#' dm_results <- random_edgeR_results()
#' peaks <- new_peaks_fn(dm_results)
#' txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#' genes <- collateGenes(genes = txdb, regions = example_regions, org.Db = org.Dm.eg.db)
#'
#' annotate_genes(peaks, genes, example_regions, max_distance = 5000)
#' # view all combinations
#' annotate_genes(peaks, genes, example_regions, max_distance = NULL)
annotate_genes <- function(peaks, genes, regions, max_distance = 5000) {
    if ((!is.data.frame(peaks) | !is.data.frame(genes) | !is.data.frame(regions)) &&
        !(inherits(peaks, "GRanges") | inherits(genes, "GRanges") | inherits(regions, "GRanges"))) {
        stop("Require data.frame of peaks, genes, and regions")
    }
    regions <- data.frame(regions) %>%
      dplyr::mutate(seqnames = paste0("chr", seqnames)) %>%
      plyranges::as_granges()
    genes_gr <- plyranges::as_granges(genes)
    peaks_gr <- plyranges::as_granges(peaks)

    combo <- ..pairGenesPeaks(genes_gr, peaks_gr)

    combo <- ..annotateStats(combo)

    combo <- ..peakRelativePosition(combo)
    combo <- combo %>%
        dplyr::mutate(num = seq_len(dplyr::n()))

    n_regions <- combo[, c("seqnames", "start", "end", "width", "num")]

    n_regions <- plyranges::find_overlaps_within(
        regions, plyranges::as_granges(n_regions)
    ) %>%
        data.frame() %>%
        dplyr::group_by(.data$num) %>%
        dplyr::summarise(n_region = dplyr::n()) %>%
        data.frame()
    combo$total_regions <- n_regions[match(combo$num, n_regions$num), "n_region"]

    combo <- combo %>%
        .[order(.$seqnames, .$peak_start, .$min_distance), ] %>%
        dplyr::group_by(.data$peak_id) %>%
        dplyr::mutate(count = seq_len(dplyr::n())) %>%
        dplyr::ungroup()

    if (!is.null(max_distance)) {
        combo <- dplyr::filter(combo, .data$min_distance <= max_distance)
    }

    closest <- ..closestGene(combo)

    combo_list <- ..annotateTop5(combo)

    cols_start <- which(colnames(combo) %in% c("seqnames", "start", "end", "width"))
    combo <- combo[, c(
        cols_start, ncol(combo),
        seq_len((utils::head(cols_start, 1) - 1)),
        (utils::tail(cols_start, 1) + 1):(ncol(combo) - 1)
    )] %>%
        .[, !(colnames(.) %in% c("from", "peak_order", "abs_TSS", "count", "num"))]

    list_results <- list(closest = closest, top_five = combo_list, all = combo)
    list_results
}
#' @export
#' @rdname annotate_genes
annotatePeaksGenes <- annotate_genes


..pairGenesPeaks <- function(genes_gr, peaks_gr) {
    pair_gene <- plyranges::pair_nearest(genes_gr, peaks_gr) %>%
        data.frame()
    pair_peak <- plyranges::pair_nearest(peaks_gr, genes_gr) %>%
        data.frame()

    pair_gene <- pair_gene %>%
        stats::setNames(c(
            "gene_seqnames", "gene_start", "gene_end", "gene_width", "gene_strand",
            "peak_seqnames", "peak_start", "peak_end", "peak_width", "peak_strand",
            colnames(.[11:ncol(.)])
        ))
    pair_peak <- pair_peak %>%
        stats::setNames(c(
            "peak_seqnames", "peak_start", "peak_end", "peak_width", "peak_strand",
            "gene_seqnames", "gene_start", "gene_end", "gene_width", "gene_strand",
            colnames(.[11:ncol(.)])
        ))

    # same order of columns
    pair_peak <- pair_peak[, colnames(pair_gene)]

    overlaps <- plyranges::pair_overlaps(genes_gr, peaks_gr) %>%
        data.frame()
    overlaps <- overlaps %>%
        stats::setNames(c(
            "gene_seqnames", "gene_start", "gene_end", "gene_width", "gene_strand",
            "peak_seqnames", "peak_start", "peak_end", "peak_width", "peak_strand",
            colnames(.[11:ncol(.)])
        ))


    combo <- rbind(pair_gene, pair_peak, overlaps) %>%
        data.frame() %>%
        dplyr::distinct(.)
    combo
}

..annotateStats <- function(combo) {
    df <- combo %>%
        dplyr::mutate(
            seqnames = .data$gene_seqnames,
            start = pmin(.data$gene_start, .data$peak_start),
            end = pmax(.data$gene_end, .data$peak_end),
            width = .data$end - .data$start + 1
        )
    df <- df %>%
        dplyr::mutate(
            peak_midpoint = (.data$peak_start + .data$peak_end) / 2,
            distance_TSS = .data$TSS - .data$peak_midpoint,
            midpoint_is = ifelse(.data$distance_TSS >= 0, "Upstream", "Downstream")
        ) %>%
        dplyr::mutate(abs_TSS = abs(.data$distance_TSS)) %>%
        dplyr::group_by(.data$peak_id) %>%
        dplyr::mutate(n_genes = dplyr::n()) %>%
        dplyr::ungroup()
    df
}

..peakRelativePosition <- function(combo) {
    df <- combo %>%
        dplyr::mutate(
            gap_ups = .data$gene_start - .data$peak_end,
            gap_dow = .data$peak_start - .data$gene_end,
            gap_st = .data$gene_start - .data$peak_start,
            gap_en = .data$peak_end - .data$gene_end
        ) %>%
        dplyr::mutate(
            position = dplyr::case_when(
                .data$gap_ups > 0 & .data$gap_ups <= .data$gap_st & .data$gap_dow < 0 & .data$gap_en < 0 ~ "Peak_upstream",
                .data$gap_ups < 0 & .data$gap_st < 0 & .data$gap_dow > 0 & .data$gap_dow <= .data$gap_en ~ "Peak_downstream",
                .data$gap_ups < 0 & .data$gap_st > 0 & .data$gap_dow < 0 & .data$gap_en < 0 ~ "Peak_overlap_upstream",
                .data$gap_ups < 0 & .data$gap_st > 0 & .data$gap_dow < 0 & .data$gap_en > 0 ~ "Peak_encompass_gene",
                .data$gap_ups < 0 & .data$gap_st < 0 & .data$gap_dow < 0 & .data$gap_en > 0 ~ "Peak_overlap_downstream",
                .data$gap_ups < 0 & .data$gap_st < 0 & .data$gap_dow < 0 & .data$gap_en < 0 ~ "Peak_within_gene"
            ),
            min_distance = dplyr::case_when(
                .data$position == "Peak_upstream" ~ .data$gap_ups,
                .data$position == "Peak_downstream" ~ .data$gap_dow,
                .data$position == "Peak_overlap_upstream" ~ .data$gap_st,
                .data$position == "Peak_overlap_downstream" ~ .data$gap_en,
                TRUE ~ 0
            )
        ) %>%
        .[, !(colnames(.) %in% c("gap_ups", "gap_dow", "gap_st", "gap_en"))]
    df
}

..closestGene <- function(combo) {
    df <- combo %>%
        dplyr::filter(.data$count == 1) %>%
        dplyr::mutate(
            gene_position = paste0(.data$gene_seqnames, ":", .data$gene_strand, ":", .data$gene_start, "-", .data$gene_end),
            peak_position = paste0(.data$peak_seqnames, ":", .data$peak_start, "-", .data$peak_end)
        )
    col_order <- c(
        "seqnames", "start", "end", "width",
        "total_regions", "n_regions_dm", "peak_id", "rank_p",
        "gene_position", "ensembl_gene_id", "gene_name", "midpoint_is",
        "position"
    )
    df <- df[, col_order]
    df
}

..annotateTop5 <- function(combo) {
    df <- combo %>%
        dplyr::filter(.data$count <= 5) %>%
        dplyr::group_by(.data$peak_id) %>%
        dplyr::mutate(
            start_region = min(.data$gene_start),
            end_region = max(.data$gene_end),
            start = pmin(.data$start_region, .data$start),
            end = pmax(.data$end_region, .data$end)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$seqnames, .data$start, .data$end, .data$peak_id, .data$rank_p, .data$n_genes) %>%
        dplyr::summarise(
            list_ensembl = toString(.data$ensembl_gene_id),
            list_gene = toString(.data$gene_name),
            position = toString(.data$position),
            distance_TSS = toString(.data$distance_TSS),
            min_distance = toString(.data$min_distance)
        ) %>%
        dplyr::ungroup()
    df
}
