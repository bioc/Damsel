##geom_peak

test_that("geom_peak: output is ggplot", {
  expect_s3_class(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")))[,c(1:6,9,7,10,8)], seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_peak.new(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")))[,c(1:6,9,7,10,8)]), regions = readRDS(test_path("fixtures", "regions.rds"))))), c("patchwork", "gg", "ggplot"))
})

test_that("geom_peak: Output is error", {
  expect_error(ggplot2::ggplot() + geom_peak.new())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_peak.new())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_peak.new(2), "data.frame of peaks is required")
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_peak.new(list(c("A", "B", 23))), "data.frame of peaks is required")
})

test_that("geom_peak: Output is message", {
  expect_message(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 430000, end_region = 450000, n_col = 1) +
                   geom_peak.new(readRDS(test_path("fixtures", "test_peaks.rds"))), "No peak data available for this region")
})

test_that("geom_peak: Output is no error", {
  expect_no_error(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")))[,c(1:6,9,7,10,8)], seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
    geom_peak.new(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")))[,c(1:6,9,7,10,8)]), regions = readRDS(test_path("fixtures", "regions.rds"))))))
  expect_no_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_peak.new(readRDS(test_path("fixtures", "test_peaks.rds"))))
  expect_no_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_peak.new(readRDS(test_path("fixtures", "test_peaks.rds")), peak.label = TRUE))
})
