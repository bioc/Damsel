##geom_gatc

test_that("geom_gatc: output is ggplot", {
  expect_s3_class(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)], seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_gatc(dplyr::mutate(regions_gatc_drosophila_dm6, seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)), c("patchwork", "gg", "ggplot"))
})

test_that("geom_gatc: Output is error", {
  expect_error(ggplot2::ggplot() + geom_gatc())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_gatc())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_gatc(c("ABC", "DEF")), "data.frame of GATC sites is required")
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_gatc(list(c("A", "B", 23))), "data.frame of GATC sites is required")
})

test_that("geom_gatc: Output is no error", {
  expect_no_error(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)], seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
    geom_gatc(dplyr::mutate(regions_gatc_drosophila_dm6, seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)))
  expect_no_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_gatc(dplyr::mutate(regions_gatc_drosophila_dm6, seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)))
})
