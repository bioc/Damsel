test_that("plot_wrap: output is error", {
    expect_error(plotWrap(
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks_new.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), start = start - 3, end = start + 4, width = end - start + 1)
    ))
    expect_error(plotWrap(
        id = "abs",
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks_new.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), start = start - 3, end = start + 4, width = end - start + 1)
    ), "Id is not in provided peaks or genes")
    expect_error(plotWrap(
        id = "ENSG12ab",
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks_new.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), start = start - 3, end = start + 4, width = end - start + 1)
    ), "Id is not in provided peaks or genes")
    expect_error(plotWrap(
        id = "8",
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks_new.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), start = start - 3, end = start + 4, width = end - start + 1)
    ))
})


test_that("plot_wrap: output is no error", {
    expect_no_error(plotWrap(
        seqnames = "chr2L", start_region = 1, end_region = 10000,
        counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")),
        dm_results.df = readRDS(test_path("fixtures", "test_results.rds")),
        peaks.df = readRDS(test_path("fixtures", "test_peaks_new.rds")),
        genes.df = readRDS(test_path("fixtures", "test_genes.rds")), txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
        gatc_sites.df = dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), start = start - 3, end = start + 4, width = end - start + 1)
    ))
})
