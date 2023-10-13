test_that("gene_ontology: Output is error", {
  expect_error(goseq_fn())
  expect_error(goseq_fn(genes = readRDS(test_path("fixtures", "test_genes.rds"))))
  expect_error(goseq_fn(peaks = list(a=c(2,3,4)), genes = readRDS(test_path("fixtures", "test_genes.rds"))))#, "Require data.frame of peaks as outputted from `aggregate_peaks")#, ignore.case = TRUE)
  expect_error(goseq_fn(peaks = readRDS(test_path("fixtures", "test_peaks.rds"))))
  expect_error(goseq_fn(peaks = readRDS(test_path("fixtures", "test_peaks.rds")), regions = 4))#, "Requires data.frame of genes as outputted from `get_biomart_genes")#, ignore.case = TRUE)
})

test_that("gene_ontology: Output is no error", {
  expect_no_error(goseq_fn(annotate_genes_new(readRDS(test_path("fixtures", "test_peaks_new.rds")), readRDS(test_path("fixtures", "test_genes.rds")), regions_gatc_drosophila_dm6)$all, readRDS(test_path("fixtures", "test_genes.rds")), regions_gatc_drosophila_dm6))
  #expect_no_error(goseq_fn(annotate_genes_new(aggregate_peaks(readRDS(test_path("fixtures", "test_results.rds")), regions_gatc_drosophila_dm6), readRDS(test_path("fixtures", "test_genes.rds")), regions = regions_gatc_drosophila_dm6)$all, readRDS(test_path("fixtures", "test_genes.rds")), regions_gatc_drosophila_dm6))
})
