##geom_genes.me

test_that("geom_genes.me: output is ggplot", {
  expect_s3_class(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)], seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                    geom_genes.me(get_biomart_genes(species = "dmelanogaster_gene_ensembl", version = 109, regions = regions_gatc_drosophila_dm6)), c("patchwork", "gg", "ggplot"))
})

test_that("geom_genes.me: Output is error", {
  expect_error(ggplot2::ggplot() + geom_genes.me())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_genes.me())
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_genes.me("AB"), "data.frame of genes is required")
  expect_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
                 geom_genes.me(matrix("AB", nrow = 5)), "data.frame of genes is required")
})

test_that("geom_genes.me: Output is no error", {
  expect_no_error(plot_counts_all_bams(process_bams(system.file("extdata", package = "Damsel"))[,c(1:6,7,10,8,11,9,12)], seqnames = "chr2R", start_region = 50000, end_region = 100000, n_col = 1) +
                    geom_genes.me(get_biomart_genes(species = "dmelanogaster_gene_ensembl", version = 109, regions = regions_gatc_drosophila_dm6)))
  expect_no_error(plot_counts_all_bams(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr3L", start_region = 7000, end_region = 9000, n_col = 1) +
                    geom_genes.me(readRDS(test_path("fixtures", "test_genes.rds"))))
})
