#!/bin/bash
for i in {0..114};do
  index=$(printf "%04d" $i)
  mkdir -p output/picked-motifs-infernal-1.1/$index
  for motif in $(ls output/picked-motifs/$index | grep '.cm$');do
    echo $index $motif
    if [ ! -s output/picked-motifs-infernal-1.1/$index/$motif ];then
       sem -j 128 "cmconvert output/picked-motifs/$index/$motif > output/picked-motifs-infernal-1.1/$index/$motif"
    fi
  done
  sem --wait
done
