#!/bin/bash
for a in BATTER RNIE transterm;do
 mkdir -p output/$a/test-set/performance
 #for d in E.C M.T B.B B.S;do
 for d in B.S;do
     cat output/$a/test-set/predictions/${d}.bed | awk 'BEGIN{FS="\t";OFS="\t";}$2>=30{p=int(($2+$3)/2);print $1,p-30,p+30,$4,$5,$6}' > output/$a/test-set/predictions/${d}.extended.bed
  for t in primary non.primary;do
   for m in RIT RDT others;do
     echo $a $d $t $m
     scripts/get-recall.py --positive dataset/term-seq/test-set/bed.u40d10/$t/$d/${m}.bed --terminator output/$a/test-set/predictions/${d}.bed -s --min-score 0 --output output/$a/test-set/performance/${d}.${t}.${m}.txt
     #scripts/join-FPR-TPR.py -t output/$a/test-set/performance/${d}.${t}.${m}.txt -f output/BATTER/test-set/predictions/${d}.FPR -o output/$a/test-set/performance/${d}.${t}.${m}.ROC
     scripts/get-recall.py --positive dataset/term-seq/test-set/bed.u40d10/$t/$d/${m}.bed --terminator output/$a/test-set/predictions/${d}.extended.bed -s --min-score 0 --output output/$a/test-set/performance/${d}.${t}.${m}.extended.txt
     #scripts/join-FPR-TPR.py -t output/$a/test-set/performance/${d}.${t}.${m}.extended.txt -f output/BATTER/test-set/predictions/${d}.FPR -o output/$a/test-set/performance/${d}.${t}.${m}.extended.ROC
   done
  done 
 done
done
