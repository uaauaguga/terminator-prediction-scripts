#!/bin/bash
for i in {0..114};do
  index=$(printf "%04d" $i)
  mkdir -p output/picked-motifs/$index
  jsub "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/get-selected-motifs.py --stockholm-directory output/candidate-sequences/by-genes-by-clades/$index --cm-directory output/candidate-motif-hits/$index --output-directory output/picked-motifs/$index -n output/hits-number/${index}.txt"
done
