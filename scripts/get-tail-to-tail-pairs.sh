#!/bin/bash
bed=$1
genome=$2
outdir=$3
mkdir -p $outdir
[ -s $outdir/plus.d.nn.bed ]  || bedtools closest -N -iu -g $genome -D ref -a <(cat $bed | awk 'BEGIN{FS="\t"}$6=="+"{print}') -b $bed  > $outdir/plus.d.nn.bed
[ -s $outdir/minus.u.nn.bed ] || bedtools closest -N -id -g $genome -D ref -a <(cat $bed | awk 'BEGIN{FS="\t"}$6=="-"{print}') -b $bed  > $outdir/minus.u.nn.bed
[ -s $outdir/converged.bed ]  || scripts/filter-gene-pairs.py -f $outdir/plus.d.nn.bed -r $outdir/minus.u.nn.bed -c $outdir/converged.bed --min-distance 48 --max-distance 512
[ -s $outdir/converged.stranded.bed ] || cat $outdir/converged.bed | awk 'BEGIN{FS="\t";OFS="\t";}{$6="+";print;$6="-";print}' > $outdir/converged.stranded.bed



