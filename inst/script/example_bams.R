# Example Drosophila DamID Bam files
#
# A subset of data from the DamID experiment in Vissers et al., (2018), GEO
# accession GSE120731. Included are the samples:
# * Dam_1: SRR7948872 which is a paired-end sample
# * Sd_2: SRR7948877 which is a single-end sample
#
# Individual samples were downloaded in fastq format from the SRA portal.
# Instructions for using `pre-fetch` to download the accessions and
# `fasterq-dump` to extract the files can be found here:
# https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
#
#
# As per Vissers et. al., (2018), the fastq files were aligned into Bam files
# using Rsubread with appropriate settings for single and paired-end files. The Bamfiles were
# sorted and indexed using Samtools.
# Alignment: Rsubread
# `buildindex(basename = "dros_ref", reference = "path/to/fasta")`
# For the single end
# `align(index = "dros_ref", readfile1 = "path/SRR7948877.fastq")`
# For the paired
# `align(index = "dros_ref", readfile1 = "path/SRR7948872_1.fastq.gz", path/SRR7948872_2.fastq.gz")`
#
# The Bam files were then sorted with `samtools sort file_in.BAM -o file_out.BAM`
#
# The sorted Bam files were filtered to retain the the region 2L:1-1000
# `samtools view -b file_out.BAM 2L:1-1000 > filt_file_out.BAM`
# before being indexed to generate a bai file with
# `samtools index filt_file_out.BAM -o filt_file_out.BAM.bai`
#
# The samples are available at:
# <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120731>
