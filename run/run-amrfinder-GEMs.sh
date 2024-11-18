#!/bin/bash
for otu_id in $(cat additional-bacteria-ids.txt );do 
  [ -s output/GEMs/${otu_id}.txt ] || jsub -R 'span[hosts=1]' -m 'ib-node153 ib-node159 ib-node161 ib-node162 ib-node164 ib-node166 ib-node170 ib-node171 ib-node173 ib-node174 ib-node176' -n 2 "export AMRFINDER_DB=/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/share/amrfinderplus/data/2023-02-23.1;export PATH=~/qhsky1/miniconda/envs/bioinfo-env/bin:$PATH;~/qhsky1/miniconda/envs/bioinfo-env/bin/amrfinder -p data/proteins/${otu_id}.faa > output/GEMs/${otu_id}.txt 2> output/GEMs/${otu_id}.log;"
done
