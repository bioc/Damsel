
#test_that("gatc_track:Output is a list", {
#  expect_s3_class(gatc_region_fn(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6), "list")
#})

test_that("gatc_track:Output is error", {
  expect_error(gatc_region_fn())
  expect_error(gatc_region_fn(list(a = c(1,2))), "Must have a BSgenome object such as BSgenome.Dmelanogaster.UCSC.dm6, OR the path to a FASTA file")#, ignore.case = TRUE)
})


test_that("gatc_track:Output is not error", {
  expect_no_error(gatc_region_fn(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6))
})
