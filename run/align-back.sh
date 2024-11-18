#!/bin/bash
mkdir -p log/back-align
for i in {51..100};do
  index=$(printf "%04d" $i)
  echo $index
  mkdir -p output/back-align/$index
  jsub -n 32 -R "span[hosts=1]" -m "ib-node153 ib-node155 ib-node156 ib-node158 ib-node159 ib-node160 ib-node161 ib-node162 ib-node163 ib-node164 ib-node165 ib-node166 ib-node170 ib-node171 ib-node173" "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/motifs-back-align.py -o output/back-align/${index} --jobs 10 -m output/picked-motifs-infernal-1.1/${index} -f output/candidate-sequences/combined.fa > log/back-align/${index}.txt 2>&1"
  done
  
