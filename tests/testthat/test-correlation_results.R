#heatmap tests
test_that("Output is a plot", {
  expect_s3_class(corr_heatmap(process_bams(system.file("extdata", package = "Damsel"))), "gg")
})

test_that("Output is error", {
  expect_error(corr_heatmap(), "data.frame of counts is required", ignore.case = TRUE)
  expect_error(corr_heatmap(df = c(1,2,3)), "data.frame of counts is required", ignore.case = TRUE)
})

test_that("Output is messgae", {
  expect_message(corr_heatmap(process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores())), "default spearman's method is used", ignore.case = TRUE)
})

test_that("Output is no error/message", {
  expect_no_error(corr_heatmap(process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores())))
  expect_no_message(corr_heatmap(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), method = "spearman"))
})

#scatterplot tests
test_that("Output is a plot", {
  expect_s3_class(corr_scatter(process_bams(system.file("extdata", package = "Damsel")), "Wing_Dam-1_S1_s.bam", "Wing_Sd-1_S2_s.bam"), "gg")
})

test_that("Output is error", {
  expect_error(corr_scatter(), "data.frame of counts is required", ignore.case = TRUE)
  expect_error(corr_scatter(df = c(1,2,3)), "data.frame of counts is required", ignore.case = TRUE)
  expect_error(corr_scatter(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), sample_1 = "Wing_Dam-1_S1_s.bam"), "sample_2 must be a character vector", ignore.case = TRUE)
  expect_error(corr_scatter(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), sample_1 = "Wing_Dam-1_S1_s.bam", sample_2 = c(1,2,3)), "sample_2 must be a character vector", ignore.case = TRUE)
  expect_error(corr_scatter(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), sample_2 = "Wing_Dam-1_S1_s.bam"), "sample_1 must be a character vector", ignore.case = TRUE)
  expect_error(corr_scatter(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), sample_1 = list(c(1,2,3)), sample_2 = "Wing_Dam-1_S1_s.bam"), "sample_1 must be a character vector", ignore.case = TRUE)
})

test_that("Output is message", {
  expect_message(corr_scatter(process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), sample_1 = "Wing_Dam-1_S1_s.bam", sample_2 = "Wing_Sd-1_S2_s.bam"), "default spearman's method is used", ignore.case = TRUE)
})

test_that("Output is no error/message", {
  expect_no_error(corr_scatter(process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), sample_1 = "Wing_Dam-1_S1_s.bam", sample_2 = "Wing_Sd-1_S2_s.bam"))
  expect_no_message(corr_scatterp(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), sample_1 = "Wing_Dam-1_S1_s.bam", sample_2 = "Wing_Sd-1_S2_s.bam", method = "spearman"))
})

