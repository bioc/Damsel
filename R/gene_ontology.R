#' Gene ontology analysis
#'
#'  `testGeneOntology` identifies the over-represented GO terms from the peak data, correcting for the number of GATC regions matching to each gene.
#'
#' @param annotation A data.frame of annotated genes and peaks as `annotate_peaks()$all`.
#' @param genes A data.frame of gene data as outputted from `get_biomart_genes()`.
#' @param regions A data.frame of GATC regions.
#' @param extend_by A number to extend the start and end of the genes. We recommend leaving to the default of 2000 bp.
#' * This is done to incorporate the acceptable distance of a peak to a gene.
#' * This also allows for consistency across significant and non-significant genes
#' @param fdr_threshold The FDR threshold used for significance in the ontology. Default is 0.05
#' @param bias Alternatively, the bias can be input by itself.
#' @export
#' @references Young MD, Wakefield MJ, Smyth GK, Oshlack A (2010). “Gene ontology analysis for RNA-seq: accounting for selection bias.” Genome Biology, 11, R14.
#' @return 3 objects
#'  * Plot of goodness of fit of model
#'  * Data frame of significant GO category results
#'  * Probability weights for each gene
#' @examples
#' library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#' library(org.Dm.eg.db)
#' set.seed(123)
#' example_regions <- random_regions()
#' peaks <- identifyPeaks(random_edgeR_results())
#'
#' txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#' genes <- collateGenes(genes = txdb, regions = example_regions, org.Db = org.Dm.eg.db)
#' annotation <- annotatePeaksGenes(peaks, genes, example_regions)$all
#'
#' ontology <- testGeneOntology(annotation, genes, example_regions)
#' ontology$signif_results
#' ontology$prob_weights
testGeneOntology <- function(annotation, genes, regions, extend_by = 2000, fdr_threshold = 0.05, bias = NULL) {
    dm_genes <- dplyr::filter(annotation, .data$min_distance <= extend_by)
    regions <- data.frame(regions)
    goseq_data <- genes
    goseq_data <- ..geneModExtend(goseq_data, regions, extend_by = {{ extend_by }})
    goseq_data <- goseq_data %>%
        dplyr::mutate(dm = ifelse(.data$ensembl_gene_id %in% dm_genes$ensembl_gene_id, 1, 0))

    gene.vector <- goseq_data$dm
    names(gene.vector) <- goseq_data$ensembl_gene_id

    if (is.null(bias)) {
        message("Bias will be n_regions that are contained within the gene length")
        bias <- goseq_data$n_regions
    }
    pwf <- goseq::nullp(gene.vector, "dm6", "ensGene", bias.data = bias)

    GO.wall <- goseq::goseq(pwf, "dm6", "ensGene")

    GO.wall <- GO.wall %>% dplyr::mutate(FDR = stats::p.adjust(.data$over_represented_pvalue, method = "BH"))
    GO.wall <- GO.wall %>%
        dplyr::filter(.data$FDR < fdr_threshold) %>%
        .[order(.$FDR, decreasing = FALSE), ]

    list <- list(signif_results = GO.wall, prob_weights = pwf)

    list
}


..geneModExtend <- function(genes, regions, extend_by=2000) {
    genes_mod <- data.frame(genes)
    genes_mod$start <- genes_mod$start - extend_by
    genes_mod$end <- genes_mod$end + extend_by
    genes_mod <- genes_mod[, !(colnames(genes_mod) %in% "width")]
    genes_mod <- genes_mod %>%
        plyranges::as_granges() %>%
        data.frame()
    regions_gr <- dplyr::mutate(regions, seqnames = paste0("chr", .data$seqnames)) %>%
        plyranges::as_granges()
    new <- plyranges::find_overlaps_within(
        regions_gr,
        plyranges::as_granges(genes_mod)
    ) %>%
        data.frame() %>%
        dplyr::group_by(.data$ensembl_gene_id) %>%
        dplyr::mutate(n_regions = dplyr::n()) %>%
        dplyr::ungroup() %>%
        data.frame()
    genes_mod$n_regions <- new[match(genes_mod$ensembl_gene_id, new$ensembl_gene_id), "n_regions"]
    genes_mod <- genes_mod %>%
        dplyr::mutate(n_regions = dplyr::coalesce(.data$n_regions, 0))
    genes_mod
}

#' Plot gene ontology results
#'
#' `plotGeneOntology()` plots the top 10 GO terms in a ggplot2 style plot.
#'
#' A dot plot with the FDR on the x-axis, the size of the dot being the number of genes in the GO category, and the colour of the dot being the ontology (Biological Process, Cellular Component, and Molecular Function).
#'
#' @param signif_results The results as outputted from goseq_fn()$signif_results. Selects the top 10 GO terms as default.
#' @param fdr_threshold The FDR threshold used for significance in the ontology. Default is 0.05
#'
#' @return A ggplot2 object
#' @export
#' @examples
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' library("org.Dm.eg.db")
#' set.seed(123)
#' example_regions <- random_regions()
#' peaks <- identifyPeaks(random_edgeR_results())
#' txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#' genes <- collateGenes(genes = txdb, regions = example_regions, org.Db = org.Dm.eg.db)
#' annotation <- annotatePeaksGenes(peaks, genes, example_regions)$all
#'
#' ontology <- testGeneOntology(annotation, genes, example_regions)$signif_results
#' plotGeneOntology(ontology)
plotGeneOntology <- function(signif_results, fdr_threshold=0.05) {
    df <- signif_results[seq_len(min(10, nrow(signif_results))), ]
    df <- df %>% dplyr::filter(!is.na(.data$category))
    max_fdr <- max(-log10(df$FDR)) + 5
    plot <- df %>%
        .[order(.$FDR, decreasing = FALSE), ] %>%
        ggplot2::ggplot(ggplot2::aes(x = -log10(.data$FDR), y = factor(.data$term, levels = .data$term), colour = .data$ontology)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$numInCat)) +
        ggplot2::geom_vline(xintercept = fdr_threshold, linetype = "dashed") +
        ggplot2::scale_x_continuous(limits = c(0, max_fdr))

    plot <- plot +
        ggplot2::labs(y = "GO category") +
        ggplot2::theme_bw()
    plot
}
