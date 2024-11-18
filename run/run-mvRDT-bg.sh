#!/bin/bash
m=PLSDA.shuffled
t=PLSDA
model=tools/mvRDT/model.PLSDA.shuffled.bg.10.ensemble.pkl

#m=gbm.shuffled
#t=gbm
#model=tools/mvRDT/model.gbm.shuffled.bg.10.ensemble.pkl

for asm_id in GCF_000005845.2 GCF_000750555.1 GCF_000006945.1 GCF_000009045.1 GCF_000195955.2 GCF_000008685.2;do

  outdir=output/mvRDT/background
  echo "prepare chunks ..."
  [ -s $outdir/${asm_id}.chunked.fa ] || scripts/prepare-chunks-for-mvRDT.py -i dataset/background/simulated/genomes/refseq-short/${asm_id}.fa -o $outdir/${asm_id}.chunked.fa

  echo "split chunks ..."
  [ -d $outdir/${asm_id}.chunks ] || scripts/split-fasta.py -i $outdir/${asm_id}.chunked.fa -o $outdir/${asm_id}.chunks

  echo "prepare features ..."

  for fasta in $(ls $outdir/${asm_id}.chunks | grep '.fa$' );do
    #[ -s $outdir/${asm_id}.chunks/${fasta%.*}/features.txt ] || jsub -m "ib-node165 ib-node166 ib-node173" "export PATH=$PATH:~/qhsky1/miniconda/envs/bioinfo-env/bin;tools/mvRDT/prepare-mvRDT-descripter.py -i $outdir/${asm_id}.chunks/$fasta -od $outdir/${asm_id}.chunks/${fasta%.*} --chunked > $outdir/${asm_id}.chunks/${fasta%.*}.log 2>&1 "
    [ ! -s $outdir/${asm_id}.chunks/${fasta%.*}/scores.${m}.txt -a -s $outdir/${asm_id}.chunks/${fasta%.*}/features.txt ] && jsub -m "ib-node165 ib-node166 ib-node173" "~/qhsky1/miniconda/envs/machine-learning/bin/python tools/mvRDT/inference.${t}.py -i $outdir/${asm_id}.chunks/${fasta%.*}/features.txt -o $outdir/${asm_id}.chunks/${fasta%.*}/scores.${m}.txt -m $model"
    #[ -s $outdir/${asm_id}.chunks/${fasta%.*}/scores.${m}.txt ] || jsub -m "ib-node158 ib-node159 ib-node162 ib-node163 ib-node166" "~/qhsky1/miniconda/envs/machine-learning/bin/python tools/mvRDT/inference.${t}.py -i $outdir/${asm_id}.chunks/${fasta%.*}/features.txt -o $outdir/${asm_id}.chunks/${fasta%.*}/scores.${m}.txt -m $model"
  done

  echo "concatenate scores ..."
  #[ -s $outdir/${asm_id}.${m}.scores.txt ] || 
  cat  $outdir/${asm_id}.chunks/*/scores.${m}.txt > $outdir/${asm_id}.${m}.scores.txt

  echo "aggregate scores ..."
  [ -s dataset/background/simulated/genomes/refseq-short/${asm_id}.fa.fai  ] || samtools faidx dataset/background/simulated/genomes/refseq-short/${asm_id}.fa
  #[ -s $outdir/${asm_id}.${m}.+.bedgraph ] || 
  tools/mvRDT/aggregate-scores.py -i $outdir/${asm_id}.${m}.scores.txt -f $outdir/${asm_id}.${m}.+.bedgraph -r $outdir/${asm_id}.${m}.-.bedgraph --length dataset/background/simulated/genomes/refseq-short/${asm_id}.fa.fai
  #[ -s $outdir/${asm_id}.${m}.bed ] || 
  tools/mvRDT/pick-windows.py -f $outdir/${asm_id}.${m}.+.bedgraph -r $outdir/${asm_id}.${m}.-.bedgraph -o $outdir/${asm_id}.${m}.bed

done
