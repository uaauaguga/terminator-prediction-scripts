#!/bin/bash
indir=dataset/background/simulated/genomes/refseq
outdir=output/transterm/background/genomes/refseq
logdir=log/transterm/background/genomes/refseq
mkdir $logdir
for fasta in $(ls $indir );do
  echo $fasta
  asm_id=$(echo $fasta | awk -F '_' '{print $1"_"$2}')
  mkdir -p $outdir/$asm_id
  [ -s $outdir/${asm_id}/transterm.max.bed ] ||  jsub -m "ib-node158 ib-node159 ib-node162 ib-node163" "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/run-transterm.py -i $indir/$fasta -o $outdir/${asm_id} > $logdir/${asm_id} 2>&1"
  if [ -s $outdir/${asm_id}/transterm.max.bed -a ! -s $outdir/${asm_id}/transterm.max.bed ];then
   scripts/pick-local-max.py -i $outdir/${asm_id}/transterm.bed -o $outdir/${asm_id}/transterm.max.bed
  fi
done
