##ggplot count_bams
test_that("plot counts all bams: output is ggplot", {
  expect_s3_class(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L", start_region = 1, end_region = 10000), "gg")
})

test_that("plot counts all bams: output is error", {
  expect_error(ggplot2::ggplot() + plot_counts_all_bams())
  expect_error(plot_counts_all_bams(7), "data.frame of counts is required")
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds"))))
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = 5), "seqnames must be element of seqnames in provided data.frame")
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L"), "numeric element for start_region is required")
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = "AB"), "numeric element for start_region is required")
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 100), "numeric element for end_region is required")
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 100, end_region = "M"), "numeric element for end_region is required")
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1000, end_region = 567), "end_region must be greater than start_region")
  expect_error(plot_counts_all_bams(df = readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2R", start_region = 100, end_region = 3000), "No data available for provided region, make the region larger")
})

test_that("plot counts all bams: output is no error", {
  expect_no_error(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L", start_region = 1, end_region = 10000))
  expect_no_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000))
})
