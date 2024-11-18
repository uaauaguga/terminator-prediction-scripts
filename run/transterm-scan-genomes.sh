#!/bin/bash
indir=dataset/genome/refseq/assemblies
outdir=output/transterm/genomes/refseq
logdir=log/output/transterm/genomes/refseq
mkdir -p $logdir $outdir
for fasta in $(ls $indir | grep -v '.fai$');do
 asm_id=$(echo $fasta | awk -F '_' '{print $1"_"$2}')
 #[ -s $outdir/${asm_id}/transterm.bed ] || jsub -m "ib-node158 ib-node159 ib-node162 ib-node163" "~/qhsky1/miniconda/envs/bioinfo-env/bin/python scripts/run-transterm.py -i $indir/$fasta -o $outdir/${asm_id} > $logdir/${asm_id}.txt 2>&1"
 [ -s $outdir/${asm_id}/transterm.bed -a ! -s $outdir/${asm_id}/transterm.max.bed ] && scripts/pick-local-max.py -i $outdir/${asm_id}/transterm.bed -o $outdir/${asm_id}/transterm.max.bed 
done
