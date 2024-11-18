#!/usr/bin/env python
import argparse
import os
import numpy as np
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('rho dependency')


def extract_score(chunk_id, otu_id):    
    best_scores = {}
    best_prediction = {}
    distances = np.zeros(500)
    bed = f"data/RUTs/{chunk_id}/{otu_id}.bed"
    if not os.path.exists(bed):
        return None, None
    with open(bed) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id = fields[0]
            p = seq_id.rfind(":")
            protein_id = seq_id[p+1:]
            seq_id = seq_id[:p]        
            fields[0] = seq_id
            start, end = int(fields[1]), int(fields[2])
            score = float(fields[6])        
            if (protein_id not in best_scores) or (best_scores[protein_id] < score):
                best_scores[protein_id] = score
                best_prediction[protein_id] = fields  
    N = 0
    with open(f"data/terminators-flanking/{chunk_id}/{otu_id}.fa.fai") as f:
        for line in f:
            N += 1
    combined_score = 0            
    for protein_id in best_prediction:
        fields = best_prediction[protein_id]
        seq_id = fields[0]
        s, e = seq_id[:-3].split(":")[-1].split("-")
        s, e = int(s), int(e)
        length = e - s
        # offset is end of the stem loop
        offset = length - 300
        start, end = int(fields[1]), int(fields[2])
        #start, end = start - offset, end - offset 
        score = float(fields[6])    
        #print(score)
        combined_score += score
        if start < 0:
            continue    
        distances[start:end] += 1  
    return combined_score/N, distances


def main():
    parser = argparse.ArgumentParser(description='extract rho dependency score')
    parser.add_argument('--chunk-id','-ci',type=str, required=True, help="chunk id to process")
    args = parser.parse_args()
    chunk_id = args.chunk_id
    genome_ids = open(f"genome-ids/{chunk_id}.txt").read().strip().split("\n")
    taxonomy = pd.read_csv("otu_taxonomy.tsv",sep="\t",index_col=0)
    fout = open(f"data/RUT-scores/{args.chunk_id}.txt","w")
    for otu_id in sorted(genome_ids,key=lambda x:int(x.split("-")[-1])):
        score, distances = extract_score(chunk_id, otu_id)     
        if score is not None:
            print(otu_id, score, taxonomy.loc[otu_id,"taxonomy"], sep="\t",file=fout)
    fout.close()
    
if __name__ == "__main__":
    main()
