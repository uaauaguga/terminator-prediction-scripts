#!/bin/bash

# iteratively group candidate sequence of each gene cluster by clade
for fasta in $(ls output/candidate-sequences/by-genes );do
  #for fasta in 0005.fa 0010.fa 0022.fa 0029.fa 0047.fa;do 
  for i in 1 2 3;do
    if [ ! -s output/candidate-sequences/by-genes${i}/$fasta ];then
      echo "processing $fasta $i ..."
      j=$((i-1))
      mkdir -p output/candidate-sequences/by-genes-by-clades${i}
      mkdir -p output/candidate-sequences/by-genes${i}
      scripts/group-candidate-sequence-by-clade.py -i output/candidate-sequences/by-genes${j}/$fasta -o output/candidate-sequences/by-genes-by-clades${i}/${fasta%%.*} -s output/candidate-sequences/by-genes${i}/$fasta > log/clade.groupping.${fasta%%.*}.${i}.log
    else
      echo "$fasta $i is already processed ."
    fi
  done
done
