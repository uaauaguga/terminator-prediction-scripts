#!/bin/bash
#for i in {51..114};do
for i in {32..40};do
  # 12 24 28 33 38 43 58 61 63 79 84 109
  index=$(printf "%04d" $i)
  mkdir -p output/candidate-motif-hits/$index log/cmsearch/$index
  jsub -n 64 -R "span[hosts=1]" -m "ib-node153 ib-node155 ib-node156 ib-node158 ib-node159 ib-node162 ib-node163 ib-node164 ib-node165 ib-node170 ib-node171 ib-node173" "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/run-cmsearch.py --jobs 128 -i output/candidate-sequences/by-genes-by-clades/$index -o output/candidate-motif-hits/$index -l log/cmsearch/$index -m output/candidate-sequences/motif-assignments/${index}.txt > log/cmsearch/${index}.log 2>&1"
done
