
test_that("Output is a data frame", {
  expect_s3_class(get_biomart_genes(species = "dmelanogaster_gene_ensembl"), "data.frame")
})

test_that("Output is a data.frame", {
  expect_s3_class(gene_annotate(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), get_biomart_genes(species = "dmelanogaster_gene_ensembl")), "data.frame")
})

test_that("Output is a data.frame", {
  expect_s3_class(gene_annotate_organised(gene_annotate(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), get_biomart_genes(species = "dmelanogaster_gene_ensembl"))), "data.frame")
})
