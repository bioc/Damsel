## edgeR set up
test_that("edgeR set up: Output is a DGEList", {
    expect_s4_class(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), "DGEList")
})

test_that("edgeR set up: Output is error", {
    expect_error(makeDGE())
    expect_error(makeDGE(counts.df = list(1, 2, 3)), "Must have data.frame of counts") # , ignore.case = TRUE)
    expect_error(makeDGE(counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")), min.cpm = "A", min.samples = 3), "min.cpm must be 1 value, recommend using default value") # , ignore.case = TRUE)
    expect_error(makeDGE(counts.df = readRDS(test_path("fixtures", "test_counts_df.rds")), min.cpm = 0.5, min.samples = c(1, 2)), "min.samples must be 1 value, recommend using default value") # , ignore.case = TRUE)
})


test_that("edgeR set up: Output is no error/message", {
    expect_no_error(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))))
})



## EdgeR results
test_that("edgeR res: Output is a data frame", {
    expect_s3_class(testDmRegions(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), regions = readRDS(test_path("fixtures", "regions.rds"))), "data.frame")
})

test_that("edgeR res: Output is error", {
    expect_error(testDmRegions())
    expect_error(testDmRegions(readRDS(test_path("fixtures", "test_dge.rds")), p.value = "A", regions = readRDS(test_path("fixtures", "regions.rds"))), "p.value must be 1 number, recommend using default value") # , ignore.case = TRUE)
    expect_error(testDmRegions(readRDS(test_path("fixtures", "test_dge.rds")), p.value = c(1, 2), regions = readRDS(test_path("fixtures", "regions.rds"))), "p.value must be 1 number, recommend using default value") # , ignore.case = TRUE)
    expect_error(testDmRegions(readRDS(test_path("fixtures", "test_dge.rds")), lfc = "%", regions = readRDS(test_path("fixtures", "regions.rds"))), "lfc must be 1 number, recommend using default value") # , ignore.case = TRUE)
    expect_error(testDmRegions(readRDS(test_path("fixtures", "test_dge.rds")), lfc = list(a = c("$", 2, "A")), regions = readRDS(test_path("fixtures", "regions.rds"))), "lfc must be 1 number, recommend using default value") # , ignore.case = TRUE)
})

test_that("edgeR res: Output is no error", {
    expect_no_error(testDmRegions(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), regions = readRDS(test_path("fixtures", "regions.rds"))))
})
