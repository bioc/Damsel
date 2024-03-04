#' Create DGE object for differential testing
#'
#' `makeDGE()` sets up the edgeR analysis for visualisation of the samples [limma::plotMDS()], and then for identifying differentially methylated regions [edgeR_results()].
#'
#' @param counts.df A data.frame generated from [process_bams]. Ensure that the samples are ordered by (Dam_1.bam, Fusion_1.bam, Dam_2.bam, Fusion_2.bam, ...).
#' @param max.width Remove large regions, default is width of 10,000. We recommend this value as the Dam can methylate GATC sites up to 5kb away from the binding site, generating a total width of 10 kb.
#' @param lib.size Library size for each sample is calculated as the sum across all rows for that sample unless otherwise specified.
#' @param min.cpm Filtering parameter, minimum counts per million (cpm) of each sample. Recommend leaving at default of 0.5.
#' @param min.samples Filtering parameter, minimum number of samples to meet the criteria of keep_a in order to retain the region in the downstream analysis. Default is 3 (assuming 6 samples).
#'
#' @return An object of class `DGEList`. Refer to [edgeR::?`DGEListClass`] for details
#' @export
#' @references Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140. doi:10.1093/bioinformatics/btp616.
#' McCarthy DJ, Chen Y, Smyth GK (2012). “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” Nucleic Acids Research, 40(10), 4288-4297. doi:10.1093/nar/gks042.
#' Chen Y, Lun ATL, Smyth GK (2016). “From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.” F1000Research, 5, 1438. doi:10.12688/f1000research.8987.2.
#' Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2024). “edgeR 4.0: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” bioRxiv. doi:10.1101/2024.01.21.576131.
#' @seealso [edgeR_results()] [process_bams()]
#' @examples
#' counts.df <- random_counts()
#'
#' makeDGE(counts.df)
makeDGE <- function(counts.df, max.width = 10000, lib.size = NULL, min.cpm = 0.5, min.samples = 3) {
    if (!is.data.frame(counts.df)) {
        stop("Must have data.frame of counts")
    }
    if (!is.numeric(min.cpm) | length(min.cpm) > 1) {
        stop("min.cpm must be 1 value, recommend using default value")
    }
    if (!is.numeric(min.samples) | length(min.samples) > 1) {
        stop("min.samples must be 1 value, recommend using default value")
    }

    counts.df <- counts.df %>% dplyr::filter(.data$width <= max.width)
    matrix <- as.matrix(counts.df[, grepl("bam", colnames(counts.df), ignore.case = TRUE)])
    rownames(matrix) <- counts.df$Position

    n_samples <- seq_len(ncol(matrix) / 2)

    group <- rep(c("Dam", "Fusion"), times = length(n_samples))

    dge <- edgeR::DGEList(matrix, lib.size = lib.size, group = group, gene = counts.df[, 2:5])

    keep <- rowSums(edgeR::cpm(dge) >= min.cpm) >= min.samples
    dge <- dge[keep, , keep.lib.sizes = FALSE]

    dge <- edgeR::calcNormFactors(dge)

    design_df <- seq_len(ncol(matrix)) %>%
        data.frame() %>%
        stats::setNames("group")
    zero_vec <- rep(0, times = length(n_samples))
    for (i in n_samples) {
        design <- replace(zero_vec, i, 1)
        design <- rep(design, each = 2)
        design_df[, ncol(design_df) + 1] <- design
    }
    design <- stats::model.matrix(~., data = design_df[, seq_len(ncol(design_df)) - 1])

    dge <- edgeR::estimateDisp(dge, robust = TRUE, design = design)

    dge
}


#' Differential testing
#'
#' `testDmRegions` calculates the differential methylation results, identifying which GATC regions have been enriched in the Fusion samples relative to the controls.
#' Refer to the following pages for further details:
#' * [edgeR::glmQLFit()]
#' * [edgeR::glmQLFTest()]
#' * [edgeR::decideTestsDGE()]
#'
#' @param dge A DGEList object as outputted from [makeDGE()].
#' @param regions A data.frame of GATC regions.
#' @param p.value A number between 0 and 1 providing the required false discovery rate (FDR). Default is 0.01.
#' @param lfc A number giving the minimum absolute log2-fold-change for significant results. Default is 1.
#' @param plot An option to plot the results using edgeR::plotSmear. Default is TRUE.
#' @return A `data.frame` of differential methylation results. Columns are: Position (chromosome-start), seqnames, start, end, width, strand, number (region number), dm (edgeR result: 0,1,NA), logFC, adjust.p, meth_status (No_signal, Upreg, Not_included).
#' If plot=TRUE, will also return a [edgeR::plotSmear()] plot of the results.
#' @export
#' @references Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140. doi:10.1093/bioinformatics/btp616.
#' McCarthy DJ, Chen Y, Smyth GK (2012). “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” Nucleic Acids Research, 40(10), 4288-4297. doi:10.1093/nar/gks042.
#' Chen Y, Lun ATL, Smyth GK (2016). “From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.” F1000Research, 5, 1438. doi:10.12688/f1000research.8987.2.
#' Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2024). “edgeR 4.0: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” bioRxiv. doi:10.1101/2024.01.21.576131.
#' @seealso [makeDGE()]
#' @examples
#' set.seed(123)
#' example_regions <- random_regions()
#' counts.df <- random_counts()
#' dge <- makeDGE(counts.df)
#'
#' dm_results <- testDmRegions(dge, regions = example_regions, p.value = 0.01, lfc = 1)
#' head(dm_results)
testDmRegions <- function(dge, regions, p.value = 0.01, lfc = 1, plot = TRUE) {
    if (!is.numeric(p.value) | length(p.value) > 1) {
        stop("p.value must be 1 number, recommend using default value")
    }
    if (!is.numeric(lfc) | length(lfc) > 1) {
        stop("lfc must be 1 number, recommend using default value")
    }
    if (missing(regions)) {
        message("GATC regions required")
    }
    group <- dge$samples$group %>% as.character()
    design <- dge$design
    fit <- edgeR::glmQLFit(dge, design = design)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
    lrt_table <- qlf$table
    lrt_table <- lrt_table %>% dplyr::mutate(
        adjust.p = stats::p.adjust(.data$PValue, method = "BH"),
        dm = ifelse(.data$logFC > lfc & .data$adjust.p < p.value, 1, 0)
    )
    if (plot == TRUE) {
        detags <- rownames(dge)[as.logical(lrt_table$dm)]
        edgeR::plotSmear(qlf, de.tags = detags, ylab = "logFC - Fusion/Dam")
    }
    lrt_table <- ..addDM(lrt_table, regions)
    lrt_table <- data.frame(lrt_table)
    lrt_table
}




#' Add dm_results to GATC regions
#'
#' Used within aggregate_peaks to add in the regions that were excluded from edgeR analysis for low counts
#'
#' @param dm_results A data.frame as outputted from [edgeR_results()]
#' @param regions A data.frame of GATC regions, default is regions_gatc_drosophila_dm6
#'
#' @return A `data.frame` of regions with added information about the dm results;
#' * dm - 1,0,NA ;
#' * logFC: 0 if dm is NA ;
#' * adjust.p: 1 if dm is NA :
#' * meth_status: Upreg, No_sig, Not_included
..addDM <- function(dm_results, regions) {
    if (!is.data.frame(dm_results)) {
        stop("Must have data frame of differential testing results from `edgeR_results")
    }
    if (!is.data.frame(regions) && !(inherits(regions, "GRanges"))) {
        stop("Regions must be provided")
    }
    results <- dm_results
    regions <- data.frame(regions)
    df <- regions %>%
        dplyr::mutate(seqnames = paste0("chr", .data$seqnames), number = seq_len(nrow(.)))
    df$dm <- results[match(df$Position, row.names(results)), "dm"]
    df$logFC <- results[match(df$Position, row.names(results)), "logFC"]
    df$PValue <- results[match(df$Position, row.names(results)), "PValue"]
    df$adjust.p <- results[match(df$Position, row.names(results)), "adjust.p"]
    df <- df %>%
        dplyr::mutate(meth_status = dplyr::case_when(
            is.na(.data$dm) ~ "Not_included",
            .data$dm == 1 ~ "Upreg",
            TRUE ~ "No_sig"
        ))
    df <- df %>%
        dplyr::mutate(
            logFC = dplyr::coalesce(.data$logFC, 0),
            PValue = dplyr::coalesce(.data$PValue, 1),
            adjust.p = dplyr::coalesce(.data$adjust.p, 1)
        ) %>%
        data.frame()
    df
}
