#!/bin/bash
odd=3
pvalue=0.001
for dataset in B.B B.S E.C M.T;do
  scripts/combine-regions.py -f output/$dataset/rho-regions.+.txt -r output/$dataset/rho-regions.-.txt -o  output/$dataset/rho-regions.bed
  cat output/$dataset/rho-regions.bed | awk -v pvalue=$pvalue -v odd=$odd 'BEGIN{FS="\t";OFS="\t";}($4>odd)&&($7<pvalue){if($4=="inf"){$4=100;}print;}' | bedtools merge -i - -s -o max,min,distinct -c 4,5,6 | bedtools slop -b 75 -i - -g genomes/fasta/${dataset}.fna.fai > output/$dataset/rho-regions.${odd}.${pvalue}.bed
  scripts/annotate-intervals.py -g genomes/bed/${dataset}.bed -b output/$dataset/rho-regions.${odd}.${pvalue}.bed -o output/$dataset/rho-regions.${odd}.${pvalue}.annotated.bed -c genomes/fasta/${dataset}.fna.fai
done
