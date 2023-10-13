#get genes
test_that("get_genes: Output is a data frame", {
  expect_s3_class(get_biomart_genes(species = "dmelanogaster_gene_ensembl", version = 109, regions = regions_gatc_drosophila_dm6), "data.frame")
})

test_that("get_genes: Output is error", {
  expect_error(get_biomart_genes())
  expect_error(get_biomart_genes(species = 2), "Species must be a character vector")#, ignore.case = TRUE)
  expect_error(get_biomart_genes(species = "dmelanogaster_gene_ensembl", regions = matrix(2, nrow=3,ncol=5)), "Regions must be a data frame")#, ignore.case = TRUE)
})

test_that("get_genes: Output is message", {
  expect_message(get_biomart_genes(species = "dmelanogaster_gene_ensembl"))
  expect_message(get_biomart_genes(species = "dmelanogaster_gene_ensembl", regions = regions_gatc_drosophila_dm6), "Default version 109 used")
  expect_message(get_biomart_genes(species = "dmelanogaster_gene_ensembl", version = 109), "Default of drosophila dm6 regions used")
})

test_that("get_genes: Output is no error", {
  expect_no_error(get_biomart_genes(species = "dmelanogaster_gene_ensembl"))
  expect_no_message(get_biomart_genes(species = "dmelanogaster_gene_ensembl", version = 109, regions = regions_gatc_drosophila_dm6))
})


#####new annotation fn

#test_that("annotate_genes.new: Output is a data.frame", {
 # expect_s3_class(annotate_genes_new(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]))), get_biomart_genes(species = "dmelanogaster_gene_ensembl"), regions = regions_gatc_drosophila_dm6), "list")
#})

test_that("annotate_genes.new: Output is error", {
  expect_error(annotate_genes())
  expect_error(annotate_genes(genes = readRDS(test_path("fixtures", "test_genes.rds"))))
  expect_error(annotate_genes(peaks = list(a=c(2,3,4)), genes = readRDS(test_path("fixtures", "test_genes.rds"))), "Require data.frame of peaks as outputted from `aggregate_peaks")#, ignore.case = TRUE)
  expect_error(annotate_genes(peaks = readRDS(test_path("fixtures", "test_peaks.rds"))))
  expect_error(annotate_genes(peaks = readRDS(test_path("fixtures", "test_peaks.rds")), genes = 4), "Requires data.frame of genes as outputted from `get_biomart_genes")#, ignore.case = TRUE)
})

test_that("annotate_genes.new: Output is no error", {
  expect_no_error(annotate_genes(readRDS(test_path("fixtures", "test_peaks_new.rds")), readRDS(test_path("fixtures", "test_genes.rds"))))
  expect_no_error(annotate_genes(aggregate_peaks(edgeR_results(edgeR_set_up(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)]))), get_biomart_genes(species = "dmelanogaster_gene_ensembl"), regions = regions_gatc_drosophila_dm6))
})
