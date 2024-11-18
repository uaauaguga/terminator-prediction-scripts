#!/bin/bash
for i in {0..114};do
  index=$(printf "%04d" $i)
  echo log/cmfinder/$index
  jsub -n 64 -R "span[hosts=1]" -m "ib-node153 ib-node155 ib-node156 ib-node158 ib-node159 ib-node160 ib-node161 ib-node162 ib-node163 ib-node164 ib-node165 ib-node166" "~/qhsky1/miniconda/bin/python scripts/run-cmfinder.py -i output/candidate-sequences/by-genes-by-clades/$index --log log/cmfinder/${index} --jobs 64 > log/cmfinder/${index}.log 2>&1"
done
