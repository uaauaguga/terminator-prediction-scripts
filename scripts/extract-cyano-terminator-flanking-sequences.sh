#!/bin/bash
mkdir -p data/cyanobacteria/terminators-flanking
for genome_id in $(cat data/cyanobacteria/genome-ids.txt);do
  echo $genome_id
  bed=data/cyanobacteria/bed.annotated/${genome_id}.bed
  genome=data/cyanobacteria/assemblies/${genome_id}.fna
  scripts/extract-terminator-flanking-sequences.py -b $bed -g $genome -f data/cyanobacteria/terminators-flanking/${genome_id}.fa
done
