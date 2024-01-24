test_that("Output is data frame", {
    expect_s3_class(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")), nthreads = 2), "data.frame")
})


test_that("Output is error", {
    expect_error(process_bams(regions = readRDS(test_path("fixtures", "regions.rds"))))
    expect_error(process_bams(path_to_bams = 2, regions = readRDS(test_path("fixtures", "regions.rds"))), "Path to bams must be a character vector") # , ignore.case = TRUE)
    expect_error(process_bams(path_to_bams = list(a = c("A", "B", "C")), regions = readRDS(test_path("fixtures", "regions.rds"))), "Path to bams must be a character vector") # , ignore.case = TRUE)
    expect_error(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = c(1, 2, 3)), "GATC region file must be a data.frame") # , ignore.case = TRUE)
    expect_error(process_bams(path_to_bams = system.file("data", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds"))), "No bam files present in path") # , ignore.case = TRUE)
})


test_that("Output is a message", {
 #throws error for next one
 expect_message(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")), nthreads = 2))
})

test_that("No error", {
    expect_no_error(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds"))))
    # expect_no_warning(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds"))))
    #expect_no_message(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = readRDS(test_path("fixtures", "regions.rds")), nthreads = 2))
})
