## geom_genes.tx

test_that("geom_genes.tx: output is ggplot", {
    expect_s3_class(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
        geom_genes.tx(collateGenes(genes = "dmelanogaster_gene_ensembl", version = 109, regions = readRDS(test_path("fixtures", "regions.rds"))), TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene), c("patchwork", "gg", "ggplot"))
})

test_that("geom_genes.tx: Output is error", {
    expect_error(ggplot2::ggplot() +
        geom_genes.tx())
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
        geom_genes.tx(txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene))
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
        geom_genes.tx("AB", TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene), "data.frame/GRanges of genes is required")
    expect_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 1, end_region = 10000, n_col = 1) +
        geom_genes.tx(matrix("AB", nrow = 5), TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene), "data.frame/GRanges of genes is required")
})

test_that("geom_genes.tx: Output is no error", {
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 50000, end_region = 100000, n_col = 1) +
        geom_genes.tx(collateGenes(genes = "dmelanogaster_gene_ensembl", version = 109, regions = readRDS(test_path("fixtures", "regions.rds"))), TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene))
    expect_no_error(plotCounts(readRDS(test_path("fixtures", "test_counts_df.rds")), seqnames = "chr2L", start_region = 7000, end_region = 9000, n_col = 1) +
        geom_genes.tx(readRDS(test_path("fixtures", "test_genes.rds")), TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene))
})
