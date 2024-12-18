#!/bin/bash
for i in {0..114};do
  index=$(printf "%04d" $i)
  [ -s output/candidate-sequences/motif-assignments/${index}.txt ] || jsub "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/filter-motifs-by-genes.py -i output/candidate-sequences/by-genes-by-clades/${index} -s output/candidate-sequences/motif-assignments/${index}.txt > log/motif-filtering/${index}.log 2>&1"
done
