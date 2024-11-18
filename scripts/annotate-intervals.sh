#!/bin/bash
chunk_id=$1
mkdir -p data/bed.annotated/$chunk_id
for genome_id in $(cat genome-ids-bacteria/${chunk_id}.txt );do
  [ -s data/bed.annotated/$chunk_id/${genome_id}.bed ] || scripts/annotate-intervals.py -g data/CDS/$chunk_id/${genome_id}.bed -b data/bed/$chunk_id/${genome_id}.bed -o data/bed.annotated/$chunk_id/${genome_id}.bed -c data/lengths/$chunk_id/${genome_id}.size 
done
