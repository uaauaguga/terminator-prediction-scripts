#!/bin/bash
indir=dataset/background/simulated/genomes/GEMs-999
outdir=output/transterm/background/genomes/GEMs 
logdir=log/transterm/background/genomes/GEMs
mkdir -p $logdir $outdir
for fasta in $(ls $indir );do
  echo $fasta
  #[ -s $outdir/${fasta%.*}/transterm.max.bed ] || jsub -m 'ib-node153 ib-node155 ib-node156 ib-node159 ib-node161 ib-node162 ib-node164 ib-node166 ib-node170 ib-node171 ib-node173 ib-node174 ib-node176' "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/run-transterm.py -i $indir/$fasta -o $outdir/${fasta%.*} > $logdir/${fasta%.*} 2>&1"
  ~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/run-transterm.py -i $indir/$fasta -o $outdir/${fasta%.*} > $logdir/${fasta%.*} 2>&1
done
