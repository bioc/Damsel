test_that("get_genes: Output is error", {
    expect_error(collateGenes())
})


test_that("get_genes: Output is no error", {
    expect_no_error(collateGenes(genes = "dmelanogaster_gene_ensembl", regions = readRDS(test_path("fixtures", "regions.rds")), version = 109))
})


##### new annotation fn

test_that("annotatePeaksGenes.new: Output is error", {
    expect_error(annotatePeaksGenes())
    expect_error(annotatePeaksGenes(genes = readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds"))))
    expect_error(annotatePeaksGenes(peaks = list(a = c(2, 3, 4)), genes = readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds"))))
    expect_error(annotatePeaksGenes(peaks = readRDS(test_path("fixtures", "test_peaks.rds")), regions = readRDS(test_path("fixtures", "regions.rds"))))
    expect_error(annotatePeaksGenes(peaks = readRDS(test_path("fixtures", "test_peaks.rds")), genes = 4, regions = readRDS(test_path("fixtures", "regions.rds"))))
})

test_that("annotatePeaksGenes.new: Output is no error", {
    expect_no_error(annotatePeaksGenes(readRDS(test_path("fixtures", "test_peaks_new.rds")), readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds"))))
    expect_no_error(annotatePeaksGenes(identifyPeaks(testDmRegions(makeDGE(readRDS(test_path("fixtures", "test_counts_df.rds"))), regions = readRDS(test_path("fixtures", "regions.rds")))), collateGenes(genes = "dmelanogaster_gene_ensembl", regions = readRDS(test_path("fixtures", "regions.rds")), version = 109), regions = readRDS(test_path("fixtures", "regions.rds"))))
})
