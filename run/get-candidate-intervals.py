#!/usr/bin/env python
from collections import defaultdict
from tqdm import tqdm
def main():
    cutoff = 0.75
    for i in range(0,10):
        path = f"output/annotation/{i}.txt"
        hit_id2group_id = {} 
        print(f"Load hits from {path} ...")
        with open(path) as f:
            for line in tqdm(f):
                fields = line.strip().split("\t")
                # OTU-3995:3300009702.a:Ga0114931_10010523_5	0010	cell division protein FtsH	0.846	11	13
                similarity = float(fields[3])
                if similarity < 0.75:
                    continue
                group_id = fields[1]
                hit_id = fields[0]
                # make sure each protein is assigned to a single cluster
                if hit_id not in hit_id2group_id:
                    hit_id2group_id[hit_id] = group_id
        print(f"Load length for {i} ...")
        lengths = {}
        with open(f"/apps/home/lulab_jinyunfan/qhsky1/metagenome-ncRNA/chunked-fasta/otus.combined.{i}.size") as f:
            for line in f:
                chrom_id, length = line.strip().split("\t")
                length = int(length)
                lengths[chrom_id] = length
        path = f"/apps/home/lulab_jinyunfan/qhsky1/metagenome-ncRNA/gene-prediction/bed.chunked/otus.combined.{i}.bed"
        fout = open(f"output/candidate-intervals/{i}.bed","w")
        print(f"get termini intervals for {i} ...")
        with open(path) as f:
            # OTU-1:3300020423.a:Ga0211525_10000157	70	649	1_1	.	+
            for line in tqdm(f):
                contig_id, start, end, name, _, strand = line.strip().split("\t")
                start, end = int(start), int(end)
                protein_id = contig_id + "_" + name.split("_")[-1]
                if protein_id in hit_id2group_id:
                    if strand == "-":
                        s = start - 128
                        e = start + 16
                    else:
                        s = end - 16
                        e = end + 128
                    s = max(0,s)
                    e = min(lengths[contig_id],e)
                    group_id = hit_id2group_id[protein_id]
                    print(contig_id, s, e, protein_id + ":" + group_id, _, strand,sep="\t",file=fout)
        fout.close()
if __name__ == "__main__":
    main()
