## geom_gatc

test_that("geom_gatc: output is ggplot", {
    expect_s3_class(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_gatc(dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)), c("patchwork", "gg", "ggplot"))
})

test_that("geom_gatc: Output is error", {
    expect_error(ggplot2::ggplot() +
        geom_gatc())
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_gatc())
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_gatc(c("ABC", "DEF")), "data.frame/GRanges object of GATC sites is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_gatc(list(c("A", "B", 23))), "data.frame/GRanges object of GATC sites is required")
})

test_that("geom_gatc: Output is no error", {
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_gatc(dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)))
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000) +
        geom_gatc(dplyr::mutate(readRDS(test_path("fixtures", "regions.rds")), seqnames = paste0("chr", seqnames), start = start - 3, end = start + 4)))
})
