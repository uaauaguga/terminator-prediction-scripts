#!/bin/bash
for i in {0..114};do
  index=$(printf "%04d" $i)
  jsub  "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/get-hits-profile.py -i output/candidate-motif-hits/$index -c output/hits-number/${index}.txt -m output/candidate-sequences/motif-assignments/${index}.txt -o output/hits-profile/${index}.txt"
done
