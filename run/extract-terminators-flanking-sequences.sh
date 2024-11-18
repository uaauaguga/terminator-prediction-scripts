#!/bin/bash
#for txt in $(ls genome-ids );do
#  chunk_id=${txt%.*}
#for chunk_id in {270..274};do
for chunk_id in 187 274;do
  echo $chunk_id
  sbatch -p Acluster -o log/extract-sequence/${chunk_id}.log -e log/extract-sequence/${chunk_id}.err --wrap="bash scripts/extract-terminator-flanking-sequences.sh $chunk_id"
done
#done

