#' Plotting genes
#'
#' @param df genes as outputted from `get_biomart_genes`
#' @param plot.space gap to next plot
#' @param plot.height height of plot
#'
#' @return ggplot layer
#' @export
#' @examples
#' path_to_bams <- system.file("extdata", package = "Damsel")
#' counts.df <- process_bams(path_to_bams,
#'                           regions = regions_gatc_drosophila_dm6,
#'                           cores = 2)
#' counts.df <- counts.df[,c(1:6,7,10,8,11,9,12)]
#' genes <- get_biomart_genes(species = "dmelanogaster_gene_ensembl",
#'                            version = 109,
#'                            regions = regions_gatc_drosophila_dm6)
#'
#' plot_counts_all_bams(counts.df,
#'                      seqnames = "chr2L",
#'                      start_region = 1,
#'                      end_region = 40000,
#'                      n_col = 1) +
#'   geom_genes.me(genes)
#'
geom_genes.me <- function(df,
                          plot.space = 0.1, plot.height = 0.3) {
  structure(list(
    df = df, plot.space = 0.1, plot.height = plot.height
  ),
  class = "genes.me"
  )
}


#' @export
ggplot_add.genes.me <- function(object, plot, object_name) {
  if(!is.data.frame(object$df)) {
    stop("data.frame of genes is required")
  }
  plot2 <- plot
  while("patchwork" %in% class(plot2)) {
    plot2 <- plot2[[1]]
  }
  plot.data <- plot2$labels$title
  plot.data <- stringr::str_split_1(plot.data, ":")
  # prepare plot range
  # the plot region are not normal, so start is minimum value
  plot.chr <- plot.data[1]
  plot.data <- stringr::str_split_1(plot.data[2], "-")
  plot.region.start <- plot.data[1] %>% as.numeric()
  plot.region.end <- plot.data[2] %>% as.numeric()

  # get parameters
  df <- object$df
  plot.space <- object$plot.space
  plot.height <- object$plot.height


  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene

  #df <- df %>% filter(seqnames == plot.chr, start >= plot.region.start,  end <= plot.region.end)
  df <- GetRegion_hack(chr = plot.chr, df = df, start = plot.region.start, end = plot.region.end)
  if(nrow(df) == 0) {
    message("No gene data available for this region")
    gene_plot <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::labs("Gene") +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::coord_cartesian(xlim = c(plot.region.start, plot.region.end)) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        #axis.line.y = ggplot2::element_blank(),
        #axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
        #axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.margin = ggplot2::margin(t = 0.1, b = 0.1)
      )
    return(patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                          gene_plot,
                          ncol = 1, heights = c(1, 0.2)))
  }
  df <- df[,c(1:3,5:ncol(df))]

  gene_plot <- ggbio::autoplot(txdb, which = plyranges::as_granges(df))@ggplot +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::coord_cartesian(xlim = c(plot.region.start, plot.region.end)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      #axis.line.y = ggplot2::element_blank(),
      #axis.text.y = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
      #axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(t = 0.1, b = 0.1)
  )
  # assemble plot
  patchwork::wrap_plots(plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
                        gene_plot,
                        ncol = 1, heights = c(1, plot.height)
  )
}
