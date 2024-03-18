#!/bin/bash

bamdir="/mnt/breast"
vcfdir="../../hawk/tumour"
refdir="../../references"

gw -b out_data/all.bwa_dodi.bam -b ${bamdir}/${1}.cram -v ${vcfdir}/${1}.vcf --track hg38_genes_protein_coding.bed ${refdir}/hg38.fa
