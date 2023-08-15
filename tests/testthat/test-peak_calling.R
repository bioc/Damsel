
test_that("Output is a data.frame", {
  expect_s3_class(add_de(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), "data.frame")
})

test_that("Output is a data.frame", {
  expect_s3_class(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), "data.frame")
})
