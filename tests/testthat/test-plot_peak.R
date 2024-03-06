## geom_peak

test_that("geom_peak: output is ggplot", {
    expect_s3_class(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
        geom_peak(identifyPeaks(testDmRegions(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), regions = readRDS(test_path("fixtures", "regions.rds"))))), c("patchwork", "gg", "ggplot"))
})

test_that("geom_peak: Output is error", {
    expect_error(ggplot2::ggplot() +
        geom_peak())
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_peak())
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_peak(2), "data.frame of peaks is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_peak(list(c("A", "B", 23))), "data.frame of peaks is required")
})

test_that("geom_peak: Output is message", {
    expect_message(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 11000, end_region = 16000) +
        geom_peak(readRDS(test_path("fixtures", "test_peaks_new.rds"))), "No data available for this region")
})

test_that("geom_peak: Output is no error", {
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_peak(identifyPeaks(testDmRegions(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), regions = readRDS(test_path("fixtures", "regions.rds"))))))
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_peak(readRDS(test_path("fixtures", "test_peaks_new.rds"))))
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_peak(readRDS(test_path("fixtures", "test_peaks_new.rds")), peak.label = TRUE))
})
