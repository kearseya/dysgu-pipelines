#!/bin/bash
bwa mem -c 1000 -A2 -B3 -O5 -E2 -T0 -L0 -D 0.25 -r 1.25 -d 200 -k 11 -a -t12 ${2} ${1} | dodi --paired False -c 1 -u 21 --ol-cost 2 --max-overlap 50000 --min-aln 50 - | samtools view -bh - | samtools sort -o ${3}; samtools index ${3}
