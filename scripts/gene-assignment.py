#!/usr/bin/env python
import re 
from collections import defaultdict
import numpy as np
import sys
from tqdm import tqdm

def main():
   
    print("Load pfam information ...") 
    pfams = set()
    cluster2pfam_ids = defaultdict(set)
    cluster2name = {}
    with open("data/pfam-considered.txt") as f:
        for line in f:
            cluster_id, gene_name, pfam = line.strip().split("\t") 
            cluster_id = cluster_id.zfill(4)
            cluster2pfam_ids[cluster_id].add(pfam)
            pfams.add(pfam)
            if cluster_id not in cluster2name:
                cluster2name[cluster_id] = gene_name

    pfam2idx = {}
    cluster2idx = {}
    for idx, pfam in enumerate(sorted(list(pfams))):
        pfam2idx[pfam] = idx
    
    idx2cluster = []
    for idx, cluster_id in enumerate(sorted(list(cluster2name.keys()))):
        cluster2idx[cluster_id] = idx
        idx2cluster.append(cluster_id)
    profiles = np.zeros((len(cluster2name),len(pfams))).astype(bool)
    
    for cluster_id in cluster2pfam_ids:
        for pfam_id in cluster2pfam_ids[cluster_id]:
            profiles[cluster2idx[cluster_id],pfam2idx[pfam_id]] = True

    print("load pfam hits ...")
    pfam_hits = defaultdict(set)
    with open(f"output/hmmsearch/{sys.argv[1]}.tbl") as f:
        for line in tqdm(f):
            if line.startswith("#"):
                continue
            line = line.strip()
            fields = re.split(r"\s+",line)
            seq_id, _, pfam_name, pfam_id, e_value = fields[:5]
            e_value = float(e_value)
            if e_value > 1:
                continue
            pfam = f"{pfam_id}-{pfam_name}"
            if pfam not in pfams:
                continue
            pfam_hits[seq_id].add(pfam)
    
    fout = open(f"output/annotation/{sys.argv[1]}.txt","w") 
    pfam_numbers = profiles.sum(axis=1)
    for seq_id in tqdm(pfam_hits):
        hits = pfam_hits[seq_id]
        profile = np.zeros((1,len(pfams))).astype(bool)
        for pfam in hits:
            profile[0,pfam2idx[pfam]] = True
        common_mask = profiles & profile      
        common = common_mask.sum(axis=1) 
        #union_mask = profiles | profile
        #union = union_mask.sum(axis=1)
        #similarity = np.round(common/union,3)
        similarity = np.round(common/pfam_numbers,3)
        top = np.argmax(similarity)
        cluster = idx2cluster[top]
        gene = cluster2name[cluster]
        print(seq_id,cluster,gene,similarity[top],common[top],pfam_numbers[top],sep="\t",file=fout)
    fout.close()
if __name__ == "__main__":
    main()
