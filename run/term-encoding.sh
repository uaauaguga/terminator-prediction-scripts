#!/bin/bash
#for fasta in $(ls dataset/terminator.seq.filtered.5);do
mkdir -p output/term.60.116
for index in RF00050-FMN RF00174-Cobalamin RF00059-TPP 0084 0069 0020 0025 0046 0012 0024 0033 ;do #  0012 0024 0033 
  echo $index
  [ -s output/term.60.116/${index}.h5 ] || scripts/encoding.py --encoder-config config/RNA.encoder.medium.json --fasta dataset/terminator.60.116/sequences95/${index}.fa --dimension 512 -m models/triplet.rfam.full.512/2.8774.pt --output output/term.60.116/${index}.h5 --device cuda:1
  #[ -s output/terminator.seq.filtered.5/${index}.h5 ] || scripts/encoding.py --encoder-config config/RNA.encoder.medium.json --fasta dataset/terminator.seq.filtered.5/${index}.fa --dimension 512 -m models/triplet.rfam.full.512/2.8774.pt --output output/terminator.seq.filtered.5/${index}.h5
done
