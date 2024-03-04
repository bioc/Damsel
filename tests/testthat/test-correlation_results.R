# heatmap tests
test_that("corr_heat: Output is a plot", {
    expect_s3_class(plotCorrHeatmap(readRDS(test_path("fixtures", "test_counts_df.rds"))), "gg")
})

test_that("corr_heat: Output is error", {
    expect_error(plotCorrHeatmap())
    expect_error(plotCorrHeatmap(df = c(1, 2, 3)), "data.frame of counts is required") # , ignore.case = TRUE)
})

test_that("corr_heat: Output is message", {
    expect_message(plotCorrHeatmap(df = readRDS(test_path("fixtures", "test_counts_df.rds"))), "default spearman's method is used") # , ignore.case = TRUE)
})

test_that("corr_heat: Output is no error/message", {
    expect_no_error(plotCorrHeatmap(readRDS(test_path("fixtures", "test_counts_df.rds"))))
    expect_no_message(plotCorrHeatmap(df = readRDS(test_path("fixtures", "test_counts_df.rds")), method = "spearman"))
})
