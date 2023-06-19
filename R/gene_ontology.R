goseq_fn <- function(regions, genes, peaks, bias=n) {
  goseq_genes_peaks <- plyranges::find_overlaps_within(plyranges::as_granges(regions), plyranges::as_granges(genes)) %>%
    as.data.frame() %>%
    group_by(ensembl_gene_id) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    mutate(de = ifelse(ensembl_gene_id %in% peaks$ensembl_gene_id, 1,0)) %>%
    distinct(ensembl_gene_id, gene_width, n, de)
  gene.vector <- goseq_genes_peaks$de
  names(gene.vector) <- goseq_genes_peaks$ensembl_gene_id

  pwf = goseq::nullp(gene.vector, "dm6", "ensGene", bias.data = goseq_genes_peaks$bias)

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


#bootstrap
goseq_fn_bootstrap <- function(regions=regions_between_gatc_dm6, genes=mutate(as.data.frame(extend_upstream(ensembl_granges, 2500)), gene_width = width), peaks) {
  goseq_genes_peaks <- find_overlaps_within(as_granges(regions), as_granges(genes)) %>% as.data.frame() %>%
    group_by(ensembl_gene_id, external_gene_name) %>%
    summarise(n_regions = n()) %>%
    mutate(de = ifelse(ensembl_gene_id %in% peaks$ensembl_gene_id, 1,0))
  gene.vector <- goseq_genes_peaks$de
  names(gene.vector) <- goseq_genes_peaks$ensembl_gene_id

  pwf = nullp(gene.vector, "dm6", "ensGene", bias.data = goseq_genes_peaks$n_regions)

  GO.wall = goseq(pwf = nullp(gene.vector, "dm6", "ensGene", bias.data = goseq_genes_peaks$n_regions), "dm6", "ensGene")

  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]

  ifelse(enriched.GO == as.character(0), "No_GO", enriched.GO)
  #need a test to see if enriched.Go is empty/doesn't have 10
  #escape clause - for nothing found in the subset
  #bootstrap 1000 times, blank df this many time

}

goseq_fn_bootstrap2 <- function(x) {
  go <- x[1:10]
  ifelse(go == as.character(0), "No_GO", go[1:10])
}


bootstrap_500_test1 <- replicate(500, goseq_fn_bootstrap2(goseq_fn_bootstrap(peaks=sample_n(peaks_harvey_wing_jan_distinct_genes, 350))))
