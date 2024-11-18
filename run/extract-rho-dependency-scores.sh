#!/bin/bash
for txt in $(ls genome-ids);do
  chunk_id=${txt%.*}
  echo $chunk_id
  scripts/extract-rho-dependency-scores.py -ci $chunk_id
done
