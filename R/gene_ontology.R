#' Gene ontology analysis
#'
#'  `goseq_fn` identifies the top 10 over-represented GO terms from the peak data, correcting for the number of GATC regions matching to each gene.
#'
#' @param annotation A data.frame of annotated genes and peaks as `annotate_peaks()$all`
#' @param genes A data.frame of gene data as outputted from `get_biomart_genes()`
#' @param regions A data.frame of GATC regions. If missing, default is `regions_gatc_drosophila_dm6`
#' @param extend_by A number to extend the start and end of the genes. We recommend leaving to the default of 2000 bp.
#' * This is done to incorporate the acceptable distance of a peak to a gene.
#' * This also allows for consistency across significant and non-significant genes
#' @param bias alternatively, the bias can be input by itself.
#' @export
#' @return 4 objects
#'  * List of top 10 over-represented GO terms across the 3 GO categories
#'  * Plot of goodness of fit of model
#'  * Plot of sample data
#'  * Plot of sample data without bias correction (should be messy)
#'  @examples
#'  set.seed(123)
#'  example_regions <- random_regions()
#'  peaks <- aggregate_peaks(random_edgeR_results())
#'  genes <- get_biomart_genes(species = "dmelanogaster_gene_ensembl",
#'                            version = 109,
#'                            regions = example_regions)
#'  annotation <- annotate_genes(peaks, genes, example_regions)$all
#'
#'  goseq_fn(annotation, genes, example_regions)
#geneOntology?
goseq_fn <- function(annotation, genes, regions=regions_gatc_drosophila_dm6, extend_by=2000, bias=NULL) {
  dm_genes <- dplyr::filter(annotation, .data$min_distance <= extend_by)
  goseq_data <- genes
  goseq_data <- gene_mod_extend(goseq_data, regions, extend_by = {{extend_by}})
  goseq_data <- goseq_data %>%
    dplyr::mutate(dm = ifelse(.data$ensembl_gene_id %in% dm_genes$ensembl_gene_id, 1, 0))

  gene.vector <- goseq_data$dm
  names(gene.vector) <- goseq_data$ensembl_gene_id

  if(is.null(bias)) {
    message("Bias will be n_regions that are contained within the gene length")
    bias <- goseq_data$n_regions
  }
  pwf <- goseq::nullp(gene.vector, "dm6", "ensGene", bias.data = bias)

  GO.wall <- goseq::goseq(pwf, "dm6", "ensGene")

  GO.samp <- goseq::goseq(pwf, "dm6", "ensGene", method="Sampling", repcnt=1000)

  plot1 <- plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
                xlab="log10(Wallenius p-values)", ylab="log10(Sampling p-values)",
                xlim=c(-3,0))
  graphics::abline(0,1,col=3,lty=2)

  GO.nobias <- goseq::goseq(pwf, "dm6", "ensGene", method="Hypergeometric")

  plot2 <- plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.wall[,1], GO.nobias[,1]), 2]),
                xlab = "log10(Wallenius p-values)", ylab = "log10(Hypergeometric p-values)",
                xlim = c(-3, 0), ylim = c(-3, 0))
  graphics::abline(0, 1, col = 3, lty = 2)

  enriched.GO <- GO.wall$category[stats::p.adjust(GO.wall$over_represented_pvalue, method = "BH")<.05]
  utils::head(enriched.GO)

  list <- list(for(go in enriched.GO[1:10]){
    print(GO.db::GOTERM[[go]])
    cat("--------------------------------------\n")
  }, utils::head(pwf), plot1, plot2)

  list
}

gene_mod_extend <- function(genes, regions, extend_by=2000) {
  genes_mod <- genes
  genes_mod$start <- genes_mod$start - extend_by
  genes_mod$end <- genes_mod$end + extend_by
  genes_mod <- genes_mod[,!(colnames(genes_mod) %in% "width")]
  genes_mod <- genes_mod %>%
    plyranges::as_granges() %>%
    data.frame()
  regions_gr <- dplyr::mutate(regions, seqnames = paste0("chr", .data$seqnames)) %>%
    plyranges::as_granges()
  new <- plyranges::find_overlaps_within(regions_gr,
                                         plyranges::as_granges(genes_mod)) %>%
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
