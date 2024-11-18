#!/bin/bash
indir=output/rhotermpred/background/refseq
outdir=output/rhotermpred/FPR/refseq
mkdir -p $outdir
for asm_id in $(ls $indir);do
  sort -k1,1 -k2,2n -o output/rhotermpred/background/refseq/$asm_id/rhoTermPred.bed output/rhotermpred/background/refseq/$asm_id/rhoTermPred.bed 
  scripts/pick-local-max.py -i output/rhotermpred/background/refseq/$asm_id/rhoTermPred.bed -o output/rhotermpred/background/refseq/$asm_id/rhoTermPred.max.bed
  #[ -s  $outdir/${asm_id}.txt ] || 
  #scripts/get-FPR.py -i $bed -o $outdir/${asm_id}.txt -rc
  scripts/get-FPR.py -i output/rhotermpred/background/refseq/$asm_id/rhoTermPred.max.bed -o  $outdir/${asm_id}.txt -rc --bins 1000
done
