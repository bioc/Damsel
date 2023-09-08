##geom_de.res

test_that("geom_de.res: output is ggplot", {
  expect_s3_class(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)], seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_de.res(add_de(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)])))), c("patchwork", "gg", "ggplot"))
})

test_that("geom_de.res: Output is error", {
  expect_error(ggplot2::ggplot() + geom_de.res())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_de.res())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_de.res(7), "data.frame of de results is required")
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_de.res(matrix("AB", nrow = 5)), "data.frame of de results is required")
})

test_that("geom_de.res: Output is no error", {
  expect_no_error(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)], seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
    geom_de.res(add_de(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)])))))
  expect_no_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_de.res(add_de(readRDS(test_path("fixtures", "test_results.rds")))))
})
