#get genes
test_that("Output is a data frame", {
  expect_s3_class(get_biomart_genes(species = "dmelanogaster_gene_ensembl"), "data.frame")
})

test_that("Output is error", {
  expect_error(get_biomart_genes(), "Species must be a character vector", ignore.case = TRUE)
  expect_error(get_biomart_genes(species = 2), "Species must be a character vector", ignore.case = TRUE)
  expect_error(get_biomart_genes(species = "dmelanogaster_gene_ensembl", regions = matrix(2, nrow=3,ncol=5)), "Regions must be a data frame", ignore.case = TRUE)
})

test_that("Output is no error", {
  expect_no_error(get_biomart_genes(species = "dmelanogaster_gene_ensembl"))
})

#annotate genes
test_that("Output is a data.frame", {
  expect_s3_class(gene_annotate(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), get_biomart_genes(species = "dmelanogaster_gene_ensembl")), "data.frame")
})

test_that("Output is error", {
  expect_error(gene_annotate(genes = get_biomart_genes(species = "dmelanogaster_gene_ensembl")), "Require data.frame of peaks as outputted from `aggregate_peaks", ignore.case = TRUE)
  expect_error(gene_annotate(peaks = list(a=c(2,3,4)), genes = get_biomart_genes(species = "dmelanogaster_gene_ensembl")), "Require data.frame of peaks as outputted from `aggregate_peaks", ignore.case = TRUE)
  expect_error(gene_annotate(peaks = aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel")))))), "Requires data.frame of genes as outputted from `get_biomart_genes", ignore.case = TRUE)
  expect_error(gene_annotate(peaks = aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), genes = 4), "Requires data.frame of genes as outputted from `get_biomart_genes", ignore.case = TRUE)
})

test_that("Output is no error", {
  expect_no_error(gene_annotate(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), get_biomart_genes(species = "dmelanogaster_gene_ensembl")))
})


#organised annotation
test_that("Output is a data.frame", {
  expect_s3_class(gene_annotate_organised(gene_annotate(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), get_biomart_genes(species = "dmelanogaster_gene_ensembl"))), "data.frame")
})

test_that("Output is error", {
  expect_error(gene_annotate_organised(), "Requires data.frame of annotated peaks as outputted from `gene_annotate", ignore.case = TRUE)
  expect_error(gene_annotate_organised(annotated_peaks = c("&", 2, "G")), "Requires data.frame of annotated peaks as outputted from `gene_annotate", ignore.case = TRUE)
})

test_that("Output is no error", {
  expect_no_error(gene_annotate_organised(gene_annotate(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))))), get_biomart_genes(species = "dmelanogaster_gene_ensembl"))))
})
