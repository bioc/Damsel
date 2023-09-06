##ggplot count_bams
test_that("plot counts all bams: output is ggplot", {
  expect_s3_class(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L", start_region = 1, end_region = 10000), "gg")
})

test_that("plot counts all bams: output is error", {
  expect_error(plot_counts_all_bams())
  expect_error(plot_counts_all_bams(7), "data.frame of counts is required")
  expect_error(plot_counts_all_bams(df = process_bams(system.file("extdata", package = "Damsel"))))
  expect_error(plot_counts_all_bams(df = process_bams(system.file("extdata", package = "Damsel")), seqnames = 5), "seqnames must be element of seqnames in provided data.frame")
  expect_error(plot_counts_all_bams(df = process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L"), "numeric element for start_region is required")
  expect_error(plot_counts_all_bams(df = process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L", start_region = "AB"), "numeric element for start_region is required")
  expect_error(plot_counts_all_bams(df = process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L", start_region = 100), "numeric element for end_region is required")
  expect_error(plot_counts_all_bams(df = process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L", start_region = 100, end_region = "M"), "numeric element for end_region is required")

})

test_that("plot counts all bams: output is no error", {
  expect_s3_class(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel")), seqnames = "chr2L", start_region = 1, end_region = 10000))
})
