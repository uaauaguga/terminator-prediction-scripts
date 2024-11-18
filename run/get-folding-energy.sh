#!/bin/bash
for species in E.C B.B B.S M.T;do
   echo $species
   fasta=genomes/fasta/${species}.fna
   bed=term-seq-sites/$species/sites.annotated.bed
   #gene=genomes/bed/B.B.bed
   scripts/evaluate-folding-energy.py -i $bed -f $fasta -o term-seq-sites/$species/sites.annotated.energy.bed 
   bed=term-seq-sites/$species/random.sites.annotated.bed
   scripts/evaluate-folding-energy.py -i $bed -f $fasta -o term-seq-sites/$species/random.sites.annotated.energy.bed
done
