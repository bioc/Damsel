
test_that("Output is data frame", {
  expect_s3_class(process_bams(system.file("extdata", package = "Damsel")), "data.frame")
})

#failed - says that expect_contains doesn't exist
#test_that("Df has expected column names", {
#  expect_contains(colnames(process_bams(system.file("extdata", package = "Damsel"))), c("Position", "seqnames", "start", "end", "width", "strand"))
#})
