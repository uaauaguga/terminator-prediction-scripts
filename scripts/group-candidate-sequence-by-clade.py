#!/usr/bin/env python

import pandas as pd
from collections import defaultdict
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='split candidate sequence by clade to reasonable size ')
    parser.add_argument('--input', '-i', type=str, required=True, help='input fasta file')
    parser.add_argument('--output', '-o', type=str, required=True, help='output directory')
    parser.add_argument('--skipped','-s', type=str,required=True, help = 'where to save skipped sequence')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    print("Load sequences ...")
    sequences = defaultdict(list)
    counts = {}
    with open(args.input) as f:
        for line in f:
            seq_id = line[1:].strip()
            sequence = next(f).strip()
            otu_id = seq_id.split(":")[0]
            sequences[otu_id].append((seq_id,sequence))


    print("Load taxonomy ...")
    table = pd.read_csv("/apps/home/lulab_jinyunfan/qhsky1/metagenome-ncRNA/otu-metadata/otu_taxonomy.tsv",sep="\t",index_col=0)
    for otu_id in sequences:
        counts[otu_id] = len(sequences[otu_id])

    min_size = 8
    max_size = 128

    table = table[table["t_domain"] == "Bacteria"]
    table = table[table.index.isin(sequences)]
   
    counts = pd.Series(counts)
    table["counts"] = counts.loc[table.index]

    table = table[~table["counts"].isna()]

    fskipped = open(args.skipped,"w")

    print("groupping sequences ...")
    for level in ['phylum', 'class', 'order','family', 'genus', 'species']:
        print(f"aggregate by {level} ...")
        if table.shape[0] == 0:
            print("no species reserved, stop .")
            break
        table = table[~table[f"t_{level}"].isna()]
        table[f"t_{level}"] = table[f"t_{level}"].map(lambda x:str(x).replace(" ","-"))
        counts_by_level = table.groupby(f"t_{level}").apply(lambda x:x["counts"].sum())
        small_ids =  counts_by_level.index[counts_by_level < min_size]
        print(f"discard {len(small_ids)} {level} which are too small .")
        for clade_id in small_ids:
            otu_ids = table[table[f"t_{level}"] == clade_id].index
            for otu_id in otu_ids:
                for seq_id, sequence in sequences[otu_id]:
                    fskipped.write(f">{seq_id}\n")
                    fskipped.write(sequence + "\n")
        table = table[~table[f"t_{level}"].isin(small_ids)]
        medium_ids = counts_by_level.index[(counts_by_level < max_size) & (counts_by_level >= min_size)]
        medium_ids = list(medium_ids) # clades with reasonable number of sequences
        to_save = table[table[f"t_{level}"].isin(medium_ids)]
        print(f"saving {len(medium_ids)} clades ...")
        for clade_id in medium_ids:
            group_id = f"{level}-{clade_id}"
            otu_ids = to_save[to_save[f"t_{level}"] == clade_id].index
            with open(f"{args.output}/{group_id}.fa","w") as fout:
                for otu_id in otu_ids:
                    for seq_id, sequence in sequences[otu_id]:
                        fout.write(f">{seq_id}\n")
                        fout.write(sequence + "\n")
        table = table[~table[f"t_{level}"].isin(medium_ids)]
    
    for otu_id in table.index:
        for seq_id, sequence in sequences[otu_id]:
            fskipped.write(f">{seq_id}\n")
            fskipped.write(sequence + "\n")
    fskipped.close()
   
   
if __name__ == "__main__":
    main()
