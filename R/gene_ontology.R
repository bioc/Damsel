#' Gene ontology analysis
#'
#'  `goseq_fn` identifies the top 10 over-represented GO terms from the peak data, correcting for the number of GATC regions matching to each gene.
#'
#' @param annotation A data.frame of annotated genes and peaks as `annotate_peaks()$all`
#' @param genes A data.frame of gene data as outputted from `get_biomart_genes()`
#' @param regions A data.frame of GATC regions.
#' @param extend_by A number to extend the start and end of the genes. We recommend leaving to the default of 2000 bp.
#' * This is done to incorporate the acceptable distance of a peak to a gene.
#' * This also allows for consistency across significant and non-significant genes
#' @param bias alternatively, the bias can be input by itself.
#' @export
#' @return 4 objects
#'  * List of top 10 over-represented GO terms across the 3 GO categories
#'  * Plot of goodness of fit of model
#'  * Data frame of significant GO category results
#'  * Probability weights for each gene
#' @examples
#' set.seed(123)
#' example_regions <- random_regions()
#' peaks <- aggregate_peaks(random_edgeR_results())
#' genes <- get_biomart_genes(
#'     species = "dmelanogaster_gene_ensembl",
#'     version = 109,
#'     regions = example_regions
#' )
#' annotation <- annotate_genes(peaks, genes, example_regions)$all
#'
#' ontology <- goseq_fn(annotation, genes, example_regions)
#' ontology$signif_results
#' ontology$prob_weights
# geneOntology?
goseq_fn <- function(annotation, genes, regions, extend_by = 2000, bias = NULL) {
    dm_genes <- dplyr::filter(annotation, .data$min_distance <= extend_by)
    goseq_data <- genes
    goseq_data <- gene_mod_extend(goseq_data, regions, extend_by = {{ extend_by }})
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
    GO.wall <- GO.wall %>% dplyr::filter(.data$FDR < 0.05)

   # go_terms <- character()
  #  for (go in GO.wall$category[seq_len(10)]) {
   #   go_terms <- c(go_terms, print(GO.db::GOTERM[[go]]))
    #}

    list <- list(#top_go = go_terms,
      signif_results=GO.wall, prob_weights = pwf)

    list
}


gene_mod_extend <- function(genes, regions, extend_by = 2000) {
    genes_mod <- genes
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
#' @param signif_results results as outputted from goseq_fn()$signif_results. Selects the top 20 GO terms as default
#' @param plot_type Plot results as a bar or a dot plot. Bar is default method
#' @param bar_x Select x axis for bar plot method= c(gene_ratio, gene_count, -log10FDR). Default is gene_ratio
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' set.seed(123)
#' example_regions <- random_regions()
#' peaks <- aggregate_peaks(random_edgeR_results())
#' genes <- get_biomart_genes(
#'     species = "dmelanogaster_gene_ensembl",
#'     version = 109,
#'     regions = example_regions
#' )
#' annotation <- annotate_genes(peaks, genes, example_regions)$all
#'
#' ontology <- goseq_fn(annotation, genes, example_regions)$signif_results
#' plot_gene_ontology(ontology, plot_type = "bar", bar_x = "gene_ratio")
#' plot_gene_ontology(ontology, plot_type = "dot")
#'
plot_gene_ontology <- function(signif_results, plot_type = c("bar", "dot"), bar_x = c("gene_ratio", "gene_count", "-log10FDR")) {
  df <- signif_results[1:20,]
  df <- df %>% dplyr::filter(!is.na(.data$category))
  if("bar" %in% plot_type) {
    if("gene_ratio" %in% bar_x) {
      plot <- df %>%
        ggplot2::ggplot(ggplot2::aes(x = numDEInCat/numInCat, y = factor(category, levels = category), fill = FDR)) +
        ggplot2::geom_bar(stat = "identity")
    } else if (bar_x == "gene_count") {
      plot <- df %>%
        ggplot2::ggplot(ggplot2::aes(x = numDEInCat, y = factor(category, levels = category), fill = FDR)) +
        ggplot2::geom_bar(stat = "identity")
    } else if(bar_x == "-log10FDR") {
      plot <- df %>%
        ggplot2::ggplot(ggplot2::aes(x = -log10(FDR), y = factor(category, levels = category), fill = FDR)) +
        ggplot2::geom_bar(stat = "identity")
    }
  } else if (plot_type == "dot") {
    plot <- df %>%
      .[order(.$numDEInCat/.$numInCat, decreasing = FALSE),] %>%
      ggplot2::ggplot(ggplot2::aes(x = numDEInCat/numInCat, y = factor(category, levels = category), colour = FDR)) +
      ggplot2::geom_point(ggplot2::aes(size = numDEInCat))
  }

  plot <- plot +
    ggplot2::labs(y = "GO category")
  plot
}
