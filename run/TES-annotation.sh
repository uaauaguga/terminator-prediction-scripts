#!/bin/bash
#for txt in $(ls genome-ids);do
#for chunk_id in  {270..274};do
#for chunk_id in 272 274;do
for chunk_id in 187;do
  #chunk_id=${txt%.*}
  sbatch --wrap="bash scripts/annotate-intervals.sh $chunk_id > data/bed.annotated/${chunk_id}.log 2>&1"
done
#done
