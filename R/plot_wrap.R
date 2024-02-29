#' Wrapper function for plotting
#'
#' `plot_wrap` plots all the available plots at once
#'
#' @param id A character vector of peak OR gene identifier(s) if wish to plot in peak/gene centric manner. Default is NULL.
#' @param seqnames A chromosome. Default is NULL.
#' @param start_region A number providing the start of region to plot. Default is NULL.
#' @param end_region A number providing the end of region to plot. Default is NULL.
#' @param counts.df A data.frame of counts as from [process_bams()]. Default is NULL.
#' @param dm_results.df A data.frame of dm results as from [edgeR_results()]. Default is NULL.
#' @param peaks.df A data.frame of peaks as from [aggregate_peaks()]. Default is NULL.
#' @param genes.df A data.frame of genes as from [get_biomart_genes()]. Default is NULL.
#' @param txdb A TxDb object as from a TxDb package. Default is NULL.
#' @param gatc_sites.df A data.frame of gatc sites as from [gatc_track()$sites]. Default is NULL.
#' @param extend_by extend region. Default is 250
#' @param ... arguments passed to geom_genes.me. Allows for adjusting of the plot appearance via gene_limits and plot.height if necessary.
#' * Default for gene_limits is NULL. If the gene is disproportionately large for the plot space, we recommend reducing the size with gene_limits = c(0,2)
#'
#' @return A `ggplot` object - or list of plots if provided multiple peaks/genes
#' @export
#'
#' @examples
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' library("org.Dm.eg.db")
#'     set.seed(123)
#'     example_regions <- random_regions()
#'     gatc_sites <- dplyr::mutate(example_regions,
#'         seqnames = paste0("chr", seqnames),
#'         start = start - 3, end = start + 4
#'     )
#'     counts.df <- random_counts()
#'     dm_results <- random_edgeR_results()
#'     peaks <- aggregate_peaks(dm_results)
#'
#'     txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#'     genes <- collateGenes(txdb, example_regions, org.Db=org.Dm.eg.db)
#'
#'     ## plot using a peak_id
#'     plot_wrap(
#'         id = peaks[1, ]$peak_id,
#'         counts.df = counts.df,
#'         dm_results.df = dm_results,
#'         peaks.df = peaks,
#'         gatc_sites.df = gatc_sites,
#'         genes.df = genes, txdb = txdb
#'     )
#'
#'     ## plot using a gene id
#'     plot_wrap(
#'         id = genes[1, ]$ensembl_gene_id,
#'         counts.df = counts.df,
#'         dm_results.df = dm_results,
#'         peaks.df = peaks,
#'         gatc_sites.df = gatc_sites,
#'         genes.df = genes, txdb = txdb
#'     )
#'
#'     ## plot providing a region
#'     plot_wrap(
#'         seqnames = "chr2L", start_region = 1, end_region = 5000,
#'         counts.df = counts.df,
#'         dm_results.df = dm_results,
#'         peaks.df = peaks,
#'         gatc_sites.df = gatc_sites,
#'         genes.df = genes, txdb = txdb
#'     )
#'
#'     ## plot multiple peaks or genes by providing a vector of id's
#'     plot_wrap(
#'         id = peaks[1:2, ]$peak_id,
#'         counts.df = counts.df,
#'         dm_results.df = dm_results,
#'         peaks.df = peaks,
#'         gatc_sites.df = gatc_sites,
#'         genes.df = genes, txdb = txdb
#'     )
#'
plot_wrap <- function(id = NULL,
                      seqnames = NULL, start_region = NULL, end_region = NULL,
                      counts.df = NULL, dm_results.df = NULL, peaks.df = NULL,
                      genes.df = NULL, txdb = NULL, gatc_sites.df = NULL, extend_by = 250, ...) {
    if (is.null(id) & is.null(seqnames) & is.null(start_region) & is.null(end_region)) {
        stop("Please provide an id (peak or ensembl_gene_id), or a region to plot (seqnames, start_region, end_region)")
    }
    if (is.null(counts.df) | is.null(dm_results.df) |
        is.null(peaks.df) | is.null(genes.df) | is.null(gatc_sites.df)) {
        stop("All data.frames must be inputted (counts.df, dm_results.df, peaks.df, genes.df, gatc_sites.df)")
    }
    if (!is.null(id)) {
        #id <- id[!is.na(id)]
        peaks.df <- data.frame(peaks.df)
        genes.df <- data.frame(genes.df)
        df <- dplyr::filter(peaks.df, .data$peak_id %in% id)
        if(nrow(df) == 0) {
            df <- dplyr::filter(genes.df, .data$ensembl_gene_id %in% id)
        }
        if(nrow(df) == 0) {
          stop("Id is not in provided peaks or genes")
        }
        chr <- as.character(df$seqnames)
        start_region <- as.numeric(df$start) - extend_by
        end_region <- as.numeric(df$end) + extend_by
    }

    if (!is.null(seqnames) && !is.null(start_region) && !is.null(end_region)) {
        chr <- as.character(seqnames)
        start_region <- as.numeric(start_region) - extend_by
        end_region <- as.numeric(end_region) + extend_by
    }

    start_region[start_region < 0] <- 0

    length_start <- length(start_region)

    list_plots <- list()
    for (i in seq_len(length_start)) {
        plot <- plot_counts_all_bams(counts.df,
            seqnames = chr[i],
            start_region = start_region[i],
            end_region = end_region[i],
            ...
        ) +
            geom_dm.res.lfc(dm_results.df) +
            geom_peak.new(peaks.df) +
            geom_gatc(gatc_sites.df) +
            geom_genes.me(genes.df, txdb, ...)

        list_plots[[i]] <- plot
    }
    list_plots
}



