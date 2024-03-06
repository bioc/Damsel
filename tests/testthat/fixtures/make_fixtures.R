##regions
#library(BSgenome.Dmelanogaster.UCSC.dm6)
#regions <- getGatcRegions(BSgenome.Dmelanogaster.UCSC.dm6)$regions
#regions <- dplyr::filter(regions, seqnames %in% c("2L", "2R", "3L", "3R", "4"))
#saveRDS(regions, "regions.rds")

##counts
#counts <- countBamInGatc(path, regions)
#saveRDS(counts, "test_counts_df.rds")

##dge
#dge <- makeDGE(counts, min.samples=2)
#saveRDS(dge, "test_dge.rds")

##dm
#dm_results <- testDmRegions(dge, regions)
#saveRDS(dm_results, test_results.rds)

##peaks
#peaks <- identifyPeaks(dm_results)
#saveRDS(peaks, "test_peaks_new.rds")

##genes
#genes <- collateGenes("dmelanogaster_gene_ensembl", regions, version=109)
#saveRDS(genes, "test_genes.rds")

##txdb
#txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#saveRDS(txdb, "txdb.rds")
