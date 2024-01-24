## new peaks fn
test_that("peaks_new:Output is a data.frame", {
    expect_s3_class(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")))[, c(1:6, 9, 7, 10, 8)]), regions = readRDS(test_path("fixtures", "regions.rds")))), "data.frame")
})

test_that("peaks_new:Output is error", {
    expect_error(aggregate_peaks())
    expect_error(aggregate_peaks(dm_results = list(a = c(1, 2), b = c("E", "F"))), "Must have data frame of differential_testing results from `edgeR_results") # , ignore.case = TRUE)
})


test_that("peaks_new:Output is not error", {
    expect_no_error(aggregate_peaks(edgeR_results(readRDS(test_path("fixtures", "test_dge.rds")), regions = readRDS(test_path("fixtures", "regions.rds")))))
})
