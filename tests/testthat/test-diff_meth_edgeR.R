
##edgeR set up
test_that("Output is a DGEList", {
  expect_s4_class(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), "DGEList")
})

test_that("Output is error", {
  expect_error(edgeR_set_up(), "Must have data.frame of counts", ignore.case = TRUE)
  expect_error(edgeR_set_up(df = list(1,2,3)), "Must have data.frame of counts", ignore.case = TRUE)
  expect_error(edgeR_set_up(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), keep_a = "A", keep_b = 3), "keep_a must be 1 value, recommend using default value", ignore.case = TRUE)
  expect_error(edgeR_set_up(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores()), keep_a = 0.5, keep_b = c(1,2)), "keep_b must be 1 value, recommend using default value", ignore.case = TRUE)
})


test_that("Output is no error/message", {
  expect_no_error(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = parallel::detectCores())))
})

##MDS plot
test_that("Output is an MDS plot", {
  expect_s4_class(edgeR_plot_mds(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), "MDS")
})

test_that("Output is error", {
  expect_error(edgeR_plot_mds(), "requires DGE as outputted from `edgeR_set_up", ignore.case = TRUE)
})

test_that("Output is no error", {
  expect_no_error(edgeR_plot_mds(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))))
})

##EdgeR results
test_that("Output is a data frame", {
  expect_s3_class(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), "data.frame")
})

test_that("Output is error", {
  expect_error(edgeR_results(), "requires DGE as outputted from `edgeR_set_up", ignore.case = TRUE)
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), p.value = "A"), "p.value must be 1 number, recommend using default value", ignore.case = TRUE)
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), p.value = c(1,2)), "p.value must be 1 number, recommend using default value", ignore.case = TRUE)
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), lfc = "%"), "lfc must be 1 number, recommend using default value", ignore.case = TRUE)
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), lfc = list(a = c("$", 2, "A"))), "lfc must be 1 number, recommend using default value", ignore.case = TRUE)
})

test_that("Output is no error", {
  expect_no_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))))
})

##results plot - fn currently doesn't work - just added contents of results fn for now
test_that("Output is an MA plot", {
  expect_s4_class(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), "MA")
})

test_that("Output is error", {
  expect_error(edgeR_results_plot(), "requires DGE as outputted from `edgeR_set_up", ignore.case = TRUE)
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), p.value = list(c(1,2,3))), "p.value must be 1 number, recommend using default value", ignore.case = TRUE)
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), p.value = "A"), "p.value must be 1 number, recommend using default value", ignore.case = TRUE)
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), lfc = "&"), "lfc must be 1 number, recommend using default value", ignore.case = TRUE)
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), lfc = data.frame(a = c("$", 2, "A"))), "lfc must be 1 number, recommend using default value", ignore.case = TRUE)
})

test_that("Output is no error", {
  expect_no_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))))
})
