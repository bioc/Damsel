## plot_wrap

# test_that("plot_wrap: output is ggplot", {
# expect_s3_class(plot_wrap(seqnames = "chr2L", start_region = 1, end_region = 10000,
#                          counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
#                         dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
#                          peaks.df = readRDS(test_path("fixtures", "test_peaks.rds")),
#                         genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#                        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)), c("list"))
# })

test_that("plot_wrap: output is error", {
    expect_error(plot_wrap(
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)
    ))
    expect_error(plot_wrap(
        id = "abs",
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)
    ), "Id is not in provided peaks or genes")
    expect_error(plot_wrap(
        id = "ENSG12ab",
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)
    ), "Id is not in provided peaks or genes")
    expect_error(plot_wrap(
        id = "8",
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)
    ))
})


test_that("plot_wrap: output is no error", {
    expect_no_error(plot_wrap(
        seqnames = "chr2L", start_region = 1, end_region = 10000,
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)
    ))
})
