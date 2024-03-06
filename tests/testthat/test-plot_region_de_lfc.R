## geom_de.res.lfc

test_that("geom_de.res.lfc: output is ggplot", {
    expect_s3_class(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_dm(testDmRegions(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), regions = readRDS(test_path("fixtures", "regions.rds")))), c("patchwork", "gg", "ggplot"))
})

test_that("geom_de.res.lfc: Output is error", {
    expect_error(ggplot2::ggplot() +
        geom_dm())
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_dm())
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_dm(7), "data.frame of dm results is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_dm(matrix("AB", nrow = 5)), "data.frame of dm results is required")
})

test_that("geom_de.res.lfc: Output is no error", {
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_dm(testDmRegions(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), regions = readRDS(test_path("fixtures", "regions.rds")))))
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_dm(readRDS(test_path("fixtures", "test_results.rds"))))
})
