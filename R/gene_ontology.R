#' Gene ontology analysis
#'
#'  `goseq_fn` identifies the top 10 over-represented GO terms from the peak data, correcting for the number of GATC regions matching to each gene.
#'
#' @param regions data frame of regions between GATC sites - load data
#' @param genes data frame of gene data
#' @param peaks data frame of annotated peak data
#'
#' @return 4 objects
#'  * List of top 10 over-represented GO terms across the 3 GO categories
#'  * Plot of goodness of fit of model
#'  * Plot of sample data
#'  * Plot of sample data without bias correction (should be messy)
#' @export
#'
#' @examples
goseq_fn <- function(regions, genes, peaks) {
  goseq_genes_peaks <- plyranges::find_overlaps_within(plyranges::as_granges(regions), plyranges::as_granges(genes)) %>%
    data.frame() %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(de = ifelse(ensembl_gene_id %in% peaks$ensembl_gene_id, 1,0)) %>%
    dplyr::distinct(ensembl_gene_id, gene_width, n, de)
  gene.vector <- goseq_genes_peaks$de
  names(gene.vector) <- goseq_genes_peaks$ensembl_gene_id

  pwf = goseq::nullp(gene.vector, "dm6", "ensGene", bias.data = goseq_genes_peaks$n)

  GO.wall = goseq::goseq(pwf, "dm6", "ensGene")

  GO.samp = goseq::goseq(pwf,"dm6","ensGene",method="Sampling",repcnt=1000)

  plot1 = plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
               xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
               xlim=c(-3,0))
  abline(0,1,col=3,lty=2)

  GO.nobias = goseq::goseq(pwf,"dm6","ensGene",method="Hypergeometric")

  plot2 = plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.wall[,1],GO.nobias[,1]),2]),
               xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)",
               xlim=c(-3,0), ylim=c(-3,0))
  abline(0,1,col=3,lty=2)

  enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]
  head(enriched.GO)

  list = list(for(go in enriched.GO[1:10]){
    print(GOTERM[[go]])
    cat("--------------------------------------\n")
  }, head(pwf), plot1, plot2)

  list
}


