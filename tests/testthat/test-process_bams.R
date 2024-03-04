test_that("Output is data frame", {
    expect_s3_class(countBamInGATC(path_to_bams = system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")), nthreads = 2), "data.frame")
})


test_that("Output is error", {
    expect_error(countBamInGATC(regions = readRDS(test_path("fixtures", "regions.rds"))))
    expect_error(countBamInGATC(path_to_bams = 2, regions = readRDS(test_path("fixtures", "regions.rds"))), "Path to bams must be a character vector") # , ignore.case = TRUE)
    expect_error(countBamInGATC(path_to_bams = list(a = c("A", "B", "C")), regions = readRDS(test_path("fixtures", "regions.rds"))), "Path to bams must be a character vector") # , ignore.case = TRUE)
    expect_error(countBamInGATC(path_to_bams = system.file("extdata", package = "Damsel"), regions = c(1, 2, 3)), "GATC region file must be a GRanges object") # , ignore.case = TRUE)
    expect_error(countBamInGATC(path_to_bams = system.file("data", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds"))), "No bam files present in path") # , ignore.case = TRUE)
})



test_that("No error", {
    expect_no_error(countBamInGATC(path_to_bams = system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds"))))
})
