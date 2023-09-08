#heatmap tests
test_that("corr_heat: Output is a plot", {
  expect_s3_class(corr_heatmap(process_bams(system.file("extdata", package = "Damsel"))), "gg")
})

test_that("corr_heat: Output is error", {
  expect_error(corr_heatmap())
  expect_error(corr_heatmap(df = c(1,2,3)), "data.frame of counts is required")#, ignore.case = TRUE)
})

test_that("corr_heat: Output is message", {
  expect_message(corr_heatmap(df = readRDS(test_path("fixtures", "test_counts_df.rds"))), "default spearman's method is used")#, ignore.case = TRUE)
})

test_that("corr_heat: Output is no error/message", {
  expect_no_error(corr_heatmap(process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = 2)))
  expect_no_message(corr_heatmap(df = readRDS(test_path("fixtures", "test_counts_df.rds")), method = "spearman"))
})

#scatterplot tests
#removed because removing fn

