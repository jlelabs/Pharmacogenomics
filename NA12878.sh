#!/bin/bash

# VCF FOR ANCESTRY, PHARMACOGENOMICS
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz

# BAM FOR PHARMACOGENOMICS
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194147/NA12878_S1.bam

# FASTA FOR PHARMACOGENOMICS
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa

# FASTA FOR HLA
wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip -d human_g1k_v37.fasta.gz

# FASTQ FOR HLA
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

# BAM FOR HLA
bwa index human_g1k_v37.fasta
bwa mem -t 40 human_g1k_v37.fasta ERR194147_1.fastq.gz ERR194147_2.fastq.gz | samtools sort -@40 -o ERR194147.bam
samtools index -@40 ERR194147.bam
