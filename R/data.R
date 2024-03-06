#' Example Drosophila DamID counts
#'
#' A subset of data from the DamID experiment in Vissers et al., (2018), GEO
#' accession GSE120731. Shown are the 2 Dam-only controls, and the 2 Scalloped
#' fusion samples. The samples have the following accessions:
#' * Dam_1: SRR7948872
#' * Sd_1: SRR7948874
#' * Dam_2: SRR7948876
#' * Sd_2: SRR7948877
#'
#' Individual samples were downloaded in fastq format from the SRA portal.
#' Instructions for using `pre-fetch` to download the accessions and
#' `fasterq-dump` to extract the files can be found here:
#' https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
#'
#' As per Vissers et. al., (2018), the fastq files were aligned into Bam files
#' using Rsubread with appropriate settings for single and paired-end files. The Bamfiles were
#' sorted and indexed using Samtools.
#' Alignment: Rsubread
#' `buildindex(basename = "dros_ref", reference = "path/to/fasta")`
#' For the single end
#' `align(index = "dros_ref", readfile1 = "path/SRR7948877.fastq")`
#' For the paired
#' `align(index = "dros_ref", readfile1 = "path/SRR7948872_1.fastq.gz", path/SRR7948872_2.fastq.gz")`
#'
#' The Bam files were then sorted with `samtools sort file_in.BAM -o file_out.BAM`
#' before being indexed with `samtools index file_out.BAM -o file_out.BAM.bai`
#'
#' The counts file was made by running `countBamInGATC()` using the above samples,
#' and a GATC region file made from:
#'  `getGatcRegions(BSgenome.Dmelanogaster.UCSC.dm6)$regions`
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
#' @usage data("dros_counts")
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120731>
"dros_counts"
