#' Wrapper function for plotting
#'
#' `plot_wrap` plots all the available plots at once
#'
#' @param peak_id peak identifier if wish to plot in peak centric manner
#' @param gene_id ensembl gene id if want to plot in gene centric manner
#' @param seqnames chromosome
#' @param start_region start of region to plot
#' @param end_region end of region to plot
#' @param counts counts.df from [process_bams()]
#' @param de_results edgeR results
#' @param peaks peaks df
#' @param genes genes df [get_biomart_genes()]
#' @param gatc_sites gatc sites
#' @param extend_by extend region, default is 250
#'
#' @return ggplot object - or list of plots if provided multiple peaks/genes
#' @export
#'
#' @examples
plot_wrap <- function(peak_id = NULL, gene_id=NULL,
                      seqnames = NULL, start_region = NULL, end_region=NULL,
                      counts = NULL, de_results = NULL, peaks = NULL,
                      genes = NULL, gatc_sites = NULL, extend_by=250) {
  if(is.null(peak_id) & is.null(gene_id) &
     is.null(seqnames) & is.null(start_region) & is.null(end_region)) {
    stop("Please provide a peak id, gene id, or a region to plot (seqnames, start_region, end_region)")
  }
  if(is.null(counts) | is.null(de_results) |
     is.null(peaks) | is.null(genes) | is.null(gatc_sites)) {
    stop("All data.frames must be inputted (counts, de_results, peaks, genes, gatc_sites)")
  }
  if(!is.null(peak_id)) {
    if(!(peak_id %in% peaks$consec_dm)) {
      stop("Peak_id is not in provided peaks data.frame")
    }
    if(length(peak_id) == 1) {
      peaks <- dplyr::filter(peaks, consec_dm == peak_id)
    }
    peaks <- dplyr::filter(peaks, consec_dm %in% peak_id)

    chr <- as.character(peaks$seqnames)
    start_region <- as.numeric(peaks$start) - extend_by
    start_region[start_region < 0] <- 0
    end_region <- as.numeric(peaks$end) + extend_by
  }
  if(!is.null(gene_id)) {
    if(!(gene_id %in% genes$ensembl_gene_id)) {
      stop("Gene_id is not in provided genes data.frame")
    }
    if(length(gene_id) == 1) {
      if(!is.null(genes)) {
        genes <- dplyr::filter(genes, ensembl_gene_id == gene_id)
      }
    }
    genes <- dplyr::filter(genes, ensembl_gene_id %in% gene_id)

    chr <- as.character(genes$seqnames)
    start_region <- as.numeric(genes$start) - extend_by
    start_region[start_region < 0] <- 0
    end_region <- as.numeric(genes$end) + extend_by
  }

  if(!is.null(seqnames) && !is.null(start_region) && !is.null(end_region)) {
    chr <- as.character(seqnames)
    start_region <- as.numeric(start_region) - extend_by
    start_region[start_region < 0] <- 0
    end_region <- as.numeric(end_region) + extend_by
  }

  length_start <- length(start_region)

  if(length_start == 1) {
    plot <- plot_counts_all_bams(counts, seqnames = chr,
                                 start_region = start_region,
                                 end_region = end_region,
                                 n_col = 1) +
      geom_de.res.lfc(de_results) +
      geom_peak.new(peaks) +
      geom_gatc(gatc_sites) +
      geom_genes.me(genes)
    return(plot)
  }
  list_plots <- list()
  for (i in 1:length(start_region)) {
    plot <- plot_counts_all_bams(counts, seqnames = chr[i],
                                 start_region = start_region[i],
                                 end_region = end_region[i],
                                 n_col = 1) +
      geom_de.res.lfc(de_results) +
      geom_peak.new(peaks) +
      geom_gatc(gatc_sites) +
      geom_genes.me(genes)

    list_plots[[i]] <- plot
  }
  list_plots

}
