#!/bin/bash
for a in BATTER RNIE transterm;do
 mkdir -p output/$a/test-set/performance
 for d in E.C M.T B.B B.S;do
  # slop predictions of all tols to 60nt, to avoid influence of systematic difference in length of the prediction
  cat output/$a/test-set/predictions/${d}.bed | awk 'BEGIN{FS="\t";OFS="\t";}$2>=30{p=int(($2+$3)/2);print $1,p-30,p+30,$4,$5,$6}' > output/$a/test-set/predictions/${d}.extended.bed
  for t in primary non.primary;do
   for m in RIT RDT others;do
     echo $a $d $t $m
     scripts/get-recall.py --positive dataset/term-seq/test-set/bed.u40d10/$t/$d/${m}.bed --terminator output/$a/test-set/predictions/${d}.bed -s --min-score 0 --output output/$a/test-set/performance/${d}.${t}.${m}.txt
   done
  done 
 done
done
