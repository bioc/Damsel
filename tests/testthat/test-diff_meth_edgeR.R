

test_that("Output is a DGEList", {
  expect_s4_class(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), "DGEList")
})

test_that("Output is an MDS plot", {
  expect_s4_class(edgeR_plot_mds(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), "MDS")
})

test_that("Output is a data frame", {
  expect_s3_class(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), "data.frame")
})

test_that("Output is an MA plot", {
  expect_s4_class(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), "MA")
})
