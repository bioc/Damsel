#add_de tests
test_that("Output is a data.frame", {
  expect_s3_class(add_de(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), "data.frame")
})

test_that("Output is error", {
  expect_error(add_de(), "Must have data frame of differential_testing results from `edgeR_results", ignore.case = TRUE)
  expect_error(add_de(de_results = matrix(3,5)), "Must have data frame of differential_testing results from `edgeR_results", ignore.case = TRUE)
  expect_error(add_de(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), regions = "A"), "Regions must be a data.frame", ignore.case = TRUE)
})

test_that("Output is not error", {
  expect_no_error(add_de(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))))
})

#peaks fn tests
test_that("Output is a data.frame", {
  expect_s3_class(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), "data.frame")
})

test_that("Output is error", {
  expect_error(aggregate_peaks(), "Must have data frame of differential_testing results from `edgeR_results", ignore.case = TRUE)
  expect_error(aggregate_peaks(dm_results = list(a = c(1,2), b = c("E", "F"))), "Must have data frame of differential_testing results from `edgeR_results", ignore.case = TRUE)
  expect_error(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), regions = 7), "Regions must be a data.frame", ignore.case = TRUE)
})

test_that("Output is not error", {
  expect_no_error(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))))
})
