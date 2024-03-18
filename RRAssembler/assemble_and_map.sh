#!/bin/bash

echo "assemble_and_map.sh: \n 1 $1 \n 2 $2 \n 3 $3 \n 4 $4 \n 5 $5"
echo $4
spades.py \
-k $4 \
-o assembly/$3 \
-s $1 \
--careful

bwa mem -a $5 assembly/$3/contigs.fasta > $2.sam
# Remove unmapped reads
samtools view -bF4 $2.sam | samtools sort -o $2.bam -
samtools index $2.bam
