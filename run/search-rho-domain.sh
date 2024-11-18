#!/bin/bash
chunk_id=$1
mkdir -p data/rho-hmm-hits/$chunk_id
for genome_id in $( cat genome-ids/${chunk_id}.txt );do
  [ -s data/rho-hmm-hits/$chunk_id/${genome_id}.tbl ] || hmmsearch --noali --tblout data/rho-hmm-hits/$chunk_id/${genome_id}.tbl reference/rho-models-NAR2024.hmm data/proteins/$chunk_id/${genome_id}.faa > /dev/null 2>&1
done
