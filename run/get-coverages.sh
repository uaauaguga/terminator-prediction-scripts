#!/bin/bash
#dataset=E.C
#asm_id=GCF_000005845.2
#strand=reverse

dataset=M.T
asm_id=GCF_000195955.2
strand=forward

#dataset=V.C
#asm_id=GCF_000006745.1
#strand=reverse

#dataset=B.B
#asm_id=GCF_000008685.2
#strand=reverse


#dataset=B.S
#asm_id=GCF_000009045.1
#strand=forward


mkdir -p output/$dataset/coverage
fai=genomes/fasta/${asm_id}.fna.fai
#indir=output/$dataset/bam.by.name
indir=output/$dataset/bam
for bam in $(ls output/$dataset/bam );do
  sample_id=${bam%.*}
  echo "processing $sample_id ..."
  [ -s genomes/fasta/${asm_id}.fna.fai ] || samtools faidx genomes/fasta/${asm_id}.fna
   scripts/get-coverage-by-bins.py -i $indir/$bam -f output/$dataset/coverage/${sample_id}.+.txt -r output/$dataset/coverage/${sample_id}.-.txt --strand $strand > output/$dataset/coverage/${sample_id}.binned.log 2>&1
done
