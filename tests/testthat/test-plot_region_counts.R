## ggplot count_bams
test_that("plot counts all bams: output is ggplot", {
    expect_s3_class(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000), "gg")
})

test_that("plot counts all bams: output is error", {
    expect_error(ggplot2::ggplot() +
        plotCounts())
    expect_error(plotCounts(7), "data.frame of counts is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds"))))
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = 5), "seqnames must be element of seqnames in provided data.frame")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L"), "numeric element for start_region is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = "AB"), "numeric element for start_region is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 100), "numeric element for end_region is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 100, end_region = "M"), "numeric element for end_region is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1000, end_region = 567), "end_region must be greater than start_region")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 50000, end_region = 70000), "No data available for provided region, make the region larger")
})

test_that("plot counts all bams: output is no error", {
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000))
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000))
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, layout = "spread"))
})
