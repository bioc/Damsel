test_that("gatc_track:Output is error", {
    expect_error(getGatcRegions())
    expect_error(getGatcRegions(list(a = c(1, 2))), "Must have a BSgenome object such as BSgenome.Dmelanogaster.UCSC.dm6, OR the path to a FASTA file") # , ignore.case = TRUE)
})


test_that("gatc_track:Output is not error", {
    expect_no_error(getGatcRegions(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6))
})
