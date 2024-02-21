.onLoad <- function(libname, pkgname) {
    # set global variables in order to avoid CHECK notes
    utils::globalVariables("regions_gatc_drosophila_dm6")
    utils::globalVariables(".")
    utils::globalVariables("meth_status")
    utils::globalVariables("dm_start")
    utils::globalVariables("dm_end")
    utils::globalVariables("n_regions_dm")
    utils::globalVariables("n_regions_not_dm")
    utils::globalVariables("dm_options")
    utils::globalVariables("logFC_match")
    utils::globalVariables("plot.region.start")
    utils::globalVariables("plot.region.end")
    utils::globalVariables("gene_limits")
    utils::globalVariables("region_pos")
    utils::globalVariables("gap_regions")
    utils::globalVariables("gap_width")
    utils::globalVariables("numDEInCat")
    utils::globalVariables("numInCat")
    utils::globalVariables("FDR")
    utils::globalVariables("category")

    invisible()
}
