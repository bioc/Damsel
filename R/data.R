#' Example Drosophila DamID counts
#'
#' A subset of data from the DamID experiment in Vissers et al., (2018), GEO
#' accession GSE120731. Shown are the 2 Dam-only controls, and the 2 Scalloped
#' fusion samples.
#'
#' Individual samples were downloaded in fastq format from the SRA portal.
#' As per Vissers et. al., (2018), the fastq files were aligned using Rsubread
#' with appropriate settings for single and paired-end files. The Bamfiles were
#' sorted and indexed using Samtools.
#'
#' @format ## `dros_counts`
#' A data frame with 383,654 rows and 10 columns:
#' \describe{
#'   \item{Position}{Chromosome and start position}
#'   \item{seqnames}{Chromosome name}
#'   \item{start, end, width}{Region information}
#'   \item{strand}{DNA strand}
#'   \item{dam_1_SRR7948872.BAM, sd_1_SRR7948874.BAM, dam_2_SRR7948876.BAM,
#'   sd_2_SRR7948877.BAM}{Sample counts}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120731>
"dros_counts"
