#!/bin/bash
while read species asm_id;do
  echo $species $asm_id
  #scripts/bedgraph2bed.py -f output/mvRDT/genomes/${asm_id}.PLSDA.shuffled.+.bedgraph -r output/mvRDT/genomes/${asm_id}.PLSDA.shuffled.-.bedgraph -o output/mvRDT/test-set/${species}.bed
  tools/mvRDT/pick-windows.py -f output/mvRDT/genomes/${asm_id}.PLSDA.shuffled.+.bedgraph -r output/mvRDT/genomes/${asm_id}.PLSDA.shuffled.-.bedgraph -o  output/mvRDT/test-set/${species}.unslopped.bed
  bedtools slop -b 50 -i output/mvRDT/test-set/${species}.unslopped.bed -g dataset/genome/refseq/test-genomes/${species}.fna.fai > output/mvRDT/test-set/${species}.bed
  cp output/mvRDT/background/${asm_id}.PLSDA.shuffled.FPR.txt output/mvRDT/test-set/${species}.FPR.txt
done < dataset/species2asm_id.txt
