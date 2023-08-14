
test_that("Output is a plot", {
  expect_s3_class(corr_heatmap(process_bams(system.file("extdata", package = "Damsel"))), "gg")
})

test_that("Output is a plot", {
  expect_s3_class(corr_scatter(process_bams(system.file("extdata", package = "Damsel")), "Dam_1", "Sd_1"), "gg")
})
