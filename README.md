
# Damsel

<!-- badges: start -->
<!-- badges: end -->

The goal of Damsel is to conduct an end-to-end analysis for DamID.
It takes an input of BAM files and a GATC region file and generates counts,
identifies methylated regions, identifies peaks, associating peaks with genes, 
conducts gene ontology testing, and provides a variety of adjustable ggplot2 
style visualisations.

## Installation

You can install the development version of Damsel like so:

``` r
BiocManager::install("Oshlack/Damsel")
```

## Example

Here we generate example results from the differential methylation results,
and plot layers of raw counts and the logFC.

``` r
library(Damsel)
set.seed(123)
example_regions <- random_regions()
example_counts <- random_counts()
example_dm <- random_edgeR_results()
head(example_dm)
plotCounts(example_counts, "chr2L", 1, 5000) +
    geom_dm(example_dm)
```


