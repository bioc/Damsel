
test_that("Output is data frame", {
  expect_s3_class(process_bams(system.file("extdata", package = "Damsel")), "data.frame")
})


