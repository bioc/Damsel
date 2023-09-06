test_that("hack_ideogram: output is ggplot", {
  expect_s3_class(hack_Ideogram(genome = "dm6", subchr = "chr2L", xlabel = TRUE, zoom.region=c(1, 20000), aspect.ratio = 0.1), "gg")
})
