
test_that("Output is data frame", {
  expect_s3_class(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, nthreads = 2), "data.frame")
})


test_that("Output is error", {
  expect_error(process_bams(regions = regions_gatc_drosophila_dm6))
  expect_error(process_bams(path_to_bams = 2, regions = regions_gatc_drosophila_dm6), "Path to bams must be a character vector")#, ignore.case = TRUE)
  expect_error(process_bams(path_to_bams = list(a=c("A", "B", "C")), regions = regions_gatc_drosophila_dm6), "Path to bams must be a character vector")#, ignore.case = TRUE)
  expect_error(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = c(1,2,3)), "GATC region file must be a data.frame")#, ignore.case = TRUE)
  expect_error(process_bams(path_to_bams = system.file("data", package = "Damsel"), regions = regions_gatc_drosophila_dm6), "No bam files present in path")#, ignore.case = TRUE)
})


test_that("Output is a message", {
  expect_message(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), nthreads = 2), "regions missing, regions_gatc_drosophila_dm6 used instead")#, ignore.case = TRUE)
  #throws error for next one
  #expect_message(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6), "cores missing, [parallel::detectCores()] used instead", ignore.case = TRUE)
})

test_that("No error/warning/message", {
  expect_no_error(process_bams(path_to_bams = system.file("extdata", package = "Damsel")))
  #expect_no_warning(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6))
  expect_no_message(process_bams(path_to_bams = system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, nthreads = 2))
})
