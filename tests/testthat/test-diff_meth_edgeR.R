
##edgeR set up
test_that("edgeR set up:Output is a DGEList", {
  expect_s4_class(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))), "DGEList")
})

test_that("edgeR set up:Output is error", {
  expect_error(edgeR_set_up())
  expect_error(edgeR_set_up(df = list(1,2,3)), "Must have data.frame of counts")#, ignore.case = TRUE)
  expect_error(edgeR_set_up(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = 1), keep_a = "A", keep_b = 3), "keep_a must be 1 value, recommend using default value")#, ignore.case = TRUE)
  expect_error(edgeR_set_up(df = process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = 1), keep_a = 0.5, keep_b = c(1,2)), "keep_b must be 1 value, recommend using default value")#, ignore.case = TRUE)
})


test_that("edgeR set up:Output is no error/message", {
  expect_no_error(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"), regions = regions_gatc_drosophila_dm6, cores = 1)))
})

##MDS plot
test_that("edgeR mds:Output is an MDS plot", {
  expect_s4_class(edgeR_plot_mds(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), "MDS")
})

test_that("edgeR mds:Output is error", {
  expect_error(edgeR_plot_mds())
})

test_that("edgeR mds:Output is no error", {
  expect_no_error(edgeR_plot_mds(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)])))
})

##EdgeR results
test_that("edgeR res:Output is a data frame", {
  expect_s3_class(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)])), "data.frame")
})

test_that("edgeR res:Output is error", {
  expect_error(edgeR_results())
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), p.value = "A"), "p.value must be 1 number, recommend using default value")#, ignore.case = TRUE)
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), p.value = c(1,2)), "p.value must be 1 number, recommend using default value")#, ignore.case = TRUE)
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), lfc = "%"), "lfc must be 1 number, recommend using default value")#, ignore.case = TRUE)
  expect_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), lfc = list(a = c("$", 2, "A"))), "lfc must be 1 number, recommend using default value")#, ignore.case = TRUE)
})

test_that("edgeR res:Output is no error", {
  expect_no_error(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)])))
})

##results plot - fn currently doesn't work - just added contents of results fn for now
#don't know why but now it's s3 and a list?? - 30/8
#test_that("edgeR res plot:Output is an MA plot", {
  #expect_s4_class(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), "MA")
  #expect_s3_class(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))), "list")
#})

test_that("edgeR res plot:Output is error", {
  expect_error(edgeR_results_plot())
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), p.value = list(c(1,2,3))), "p.value must be 1 number, recommend using default value")#, ignore.case = TRUE)
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), p.value = "A"), "p.value must be 1 number, recommend using default value")#, ignore.case = TRUE)
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), lfc = "&"), "lfc must be 1 number, recommend using default value")#, ignore.case = TRUE)
  expect_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]), lfc = data.frame(a = c("$", 2, "A"))), "lfc must be 1 number, recommend using default value")#, ignore.case = TRUE)
})

test_that("edgeR res plot:Output is no error", {
  expect_no_error(edgeR_results_plot(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)])))
})
