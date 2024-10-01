# plot counts distrib

test_that("countsDistrib: output is ggplot", {
  expect_s3_class(plotCountsDistribution(readRDS(test_path("fixtures", "test_counts_df.rds"))), c("gg", "ggplot"))
})

test_that("countsDistrib: output is error", {
  expect_error(plotCountsDistribution(data.frame(c = c(1,2), b = c("A", "B"))))
})

test_that("countsDistrib: output is no error", {
  expect_no_error(plotCountsDistribution(readRDS(test_path("fixtures", "test_counts_df.rds"))))
})


# peaks in reads

test_that("readsPeaks: output is ggplot", {
  expect_s3_class(plotCountsInPeaks(readRDS(test_path("fixtures", "test_counts_df.rds")),
                                    readRDS(test_path("fixtures", "test_results.rds")),
                                    readRDS(test_path("fixtures", "test_peaks_new.rds"))), c("gg", "ggplot"))
})

test_that("readsPeaks: output is error", {
  expect_error(plotCountsInPeaks(data.frame(c = c(1,2), b = c("A", "B"))))
})

test_that("readsPeaks: output is no error", {
  expect_no_error(plotCountsInPeaks(readRDS(test_path("fixtures", "test_counts_df.rds")),
                                    readRDS(test_path("fixtures", "test_results.rds")),
                                    readRDS(test_path("fixtures", "test_peaks_new.rds"))))
})


