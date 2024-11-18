#!/bin/bash
#indir=dataset/term-seq/test-set/bedU32D8
indir=dataset/term-seq/test-set/bedU100D100.downstream
outdir=output/rhotermpred/performance/test-set/U100D100.downstream
mkdir -p $outdir
for t in others RDT;do
  mkdir -p $outdir/$t
  for bed in $(ls $indir/$t );do
    asm_id=${bed%.*}
    path=output/rhotermpred/prediction/refseq/$asm_id/rhoTermPred.bed
    scripts/get-recall.py --positive $indir/$t/$bed --terminator $path --output $outdir/$t/${bed%.*}.txt -ms 0
    scripts/join-FPR-TPR.py --fpr output/rhotermpred/FPR/refseq/${asm_id}.txt  --tpr $outdir/$t/${bed%.*}.txt -o $outdir/$t/${asm_id}.ROC.txt 
    echo "#"
    echo "#"
done
done
echo "see result in $outdir"
