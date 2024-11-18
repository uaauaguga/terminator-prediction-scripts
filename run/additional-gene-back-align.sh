#!/bin/bash
indir=output/back-align-motif-hit-counts-by-otu-top32-chunked
outdir=output/sampled-GEM-genes/back-align
mkdir -p output/sampled-GEM-genes/back-align
for txt in $(ls $indir);do
   jsub -n 4 -R "span[hosts=1]" "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/motifs-additional-back-align.py --motif-to-sequence $indir/$txt -o output/sampled-GEM-genes/back-align/${txt%%.*} -j 16 > log/additional-back-align/$txt 2>&1"
done
