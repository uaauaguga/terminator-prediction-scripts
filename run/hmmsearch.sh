#!/bin/bash
for i in 6;do #{0..9};do
  echo $i
  jsub -n 64 "~/qhsky1/miniconda/envs/bioinfo-env/bin/hmmsearch --cpu 128 --noali --tblout output/hmmsearch/${i}.tbl -A output/hmmsearch/${i}.stk data/term.seq.site.3p.hmm /apps/home/lulab_jinyunfan/qhsky1/metagenome-ncRNA/gene-prediction/proteins.chunked/otus.combined.${i}.faa > output/hmmsearch/${i}.txt 2>&1"
done
