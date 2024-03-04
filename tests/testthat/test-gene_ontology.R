test_that("gene_ontology: Output is error", {
    expect_error(testGeneOntology())
    expect_error(testGeneOntology(genes = readRDS(test_path("fixtures", "test_genes.rds"))))
    expect_error(testGeneOntology(peaks = list(a = c(2, 3, 4)), genes = readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds")))) # , "Require data.frame of peaks as outputted from `aggregate_peaks")#, ignore.case = TRUE)
    expect_error(testGeneOntology(peaks = readRDS(test_path("fixtures", "test_peaks.rds"))))
    expect_error(testGeneOntology(peaks = readRDS(test_path("fixtures", "test_peaks.rds")), regions = 4)) # , "Requires data.frame of genes as outputted from `get_biomart_genes")#, ignore.case = TRUE)
})

test_that("gene_ontology: Output is no error", {
    expect_no_error(testGeneOntology(annotatePeaksGenes(readRDS(test_path("fixtures", "test_peaks_new.rds")), readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds")))$all, readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds"))))
})

## plot

test_that("gene_ontology_plot: Output is no error", {
    expect_no_error(plotGeneOntology(testGeneOntology(annotatePeaksGenes(readRDS(test_path("fixtures", "test_peaks_new.rds")), readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds")))$all, readRDS(test_path("fixtures", "test_genes.rds")), regions = readRDS(test_path("fixtures", "regions.rds")))$signif_results))
})
