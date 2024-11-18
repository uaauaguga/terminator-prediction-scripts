#!/bin/bash
mkdir -p log/calibration
for i in {51..100};do
  index=$(printf "%04d" $i)
  mkdir -p log/calibration/$index
  for motif in $(ls output/picked-motifs-infernal-1.1/$index );do
     [ -f log/calibration/$index/${motif} ] || jsub -R "span[hosts=1]" -m "ib-node153 ib-node155 ib-node156 ib-node158 ib-node159 ib-node160 ib-node161 ib-node162 ib-node163 ib-node164 ib-node165 ib-node166 ib-node170" -n 2 "~/qhsky1/miniconda/envs/bioinfo-env/bin/cmcalibrate --cpu 4 output/picked-motifs-infernal-1.1/$index/$motif && touch log/calibration/$index/${motif}"
  done
done
