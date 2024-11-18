#!/bin/bash
indir=dataset/genome/refseq/assemblies-short
outdir=output/rhotermpred/prediction/refseq
mkdir -p $outdir
for fasta in $(ls $indir | grep -v '.fai$' );do
  asm_id=${fasta%.*}
  [ -s output/rhotermpred/prediction/refseq/$asm_id/rhoTermPred.bed ] || jsub -m 'ib-node158 ib-node162 ib-node164 ib-node173' -n 1 "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/RhoTermPredict_algorithm-tidy.py -i $indir/${asm_id}.fna -o $outdir/$asm_id"
done
