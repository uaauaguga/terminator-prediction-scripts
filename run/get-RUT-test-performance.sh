#!/bin/bash
for a in BS-RUT  mvRDT rhotermpred;do
 mkdir -p output/$a/test-set/performance
 #for d in E.C M.T B.B B.S;do
  for d in B.S;do
  cat output/$a/test-set/predictions/${d}.bed | awk 'BEGIN{FS="\t";OFS="\t";}$2>=100{p=int(($2+$3)/2);print $1,p-100,p+100,$4,$5,$6}' > output/$a/test-set/predictions/${d}.extended.200.bed
  for t in primary non.primary;do
   for m in RIT RDT others;do
     scripts/get-recall.py --positive dataset/term-seq/test-set/bed.u40d10/$t/$d/${m}.bed --terminator output/$a/test-set/predictions/${d}.extended.200.bed -s --min-score 0 --output output/$a/test-set/performance/${d}.${t}.${m}.extended.200.txt
   done
  done 
 done
done
