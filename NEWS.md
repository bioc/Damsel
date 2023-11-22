# Damsel 0.4.0
## New features
* process_bams() new count method using `Rsubread::featureCounts()` allowing for fractional counts and differentiating between single and paired end BAM files
* edgeR_set_up() filter out large regions (> 10kb)
* aggregate_peaks() new method for ranking peaks based on theory of `csaw::getBestTest()`
* aggregate_peaks() retain small peaks so that they are able to be combined with the gaps fn
* gatc_track() simplified for one argument input - identifies if input is BSgenome or FASTA

# Damsel 0.3.1
## New features
* Plotting options - counts layout and log2 scale, and peak_id text
* New default for differential testing: p-value set to 0.01

# Damsel 0.3.0
## Bug fixes
* streamlined outputs
* removed obselete functions

# Damsel 0.2.0
## New features
* New peak output and functionality (combines peaks with small gaps)
* New genes output - with list of 3 data.frames
* Gene ontology fn 
* Plot_wrap() - plot all at once
* gatc_track() - create the GATC track

## Bug fixes
* accurate plotting

## Documentation
* runnable examples

# Damsel 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added geom_genes.me allowing for genes plot
* Fixed geom_peak.new
