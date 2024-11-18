#!/bin/bash
#m=gbm
m=PLSDA
outdir=output/mvRDT/background
for asm_id in GCF_000005845.2 GCF_000006945.1 GCF_000008685.2 GCF_000009045.1 GCF_000195955.2 GCF_000750555.1;do
  #bed=$outdir/${asm_id}.PLSDA.shuffled.bed
  bed=$outdir/${asm_id}.${m}.shuffled.bed
  #[ -s  $outdir/${asm_id}.FPR.txt ] || 
  scripts/get-FPR.py -i $bed -o $outdir/${asm_id}.${m}.shuffled.FPR.txt -rc
done
