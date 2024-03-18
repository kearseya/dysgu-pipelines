#!/bin/bash

for i in ${1}/; 
do 
	echo ${i};
	pdfunite ${i}*.pdf ${1}/$(basename ${i} .vcf).pdf; 
done
