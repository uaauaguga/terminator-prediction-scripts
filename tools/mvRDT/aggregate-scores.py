#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='combined scores by 50nt bins')
    parser.add_argument('--input', '-i', type=str, required=True, help='features')
    parser.add_argument('--forward', '-f', type=str, required=True, help='forward scores')
    parser.add_argument('--reverse', '-r', type=str, required=True, help='reverse scores')
    parser.add_argument('--length', '-l', type=str, required=True, help='length')
    args = parser.parse_args()

    ffwd = open(args.forward, "w")
    frev = open(args.reverse, "w")
    fwd_scores, rev_scores = defaultdict(float), defaultdict(float)
    fwd_counts, rev_counts = defaultdict(int), defaultdict(int)

    lengths = {}
    with open(args.length) as f:
        for line in f:
            seq_id, length = line.strip().split("\t")[:2]
            length = int(length)
            lengths[seq_id] = length

    entries = []
    with open(args.input) as f:
        for line in f:
            #Seq34_00005:CCGG_801_1300_R     0.0
            seq_id, score = line.strip().split("\t")
            score = float(score)
            p = seq_id.find("_")
            seq_id = seq_id[p+1:]
            fields = seq_id.split("_")
            start, end, direction = fields[-3:]
            start, end = int(start)-1, int(end)
            binidx = int(start/50)
            seq_id = "_".join(fields[:-3])
            if direction == "F":
                # forward direction
                for i in range(10):
                    if 50*(binidx+i) >= lengths[seq_id]:
                        continue
                    fwd_scores[(seq_id, binidx + i)] += score
                    fwd_counts[(seq_id, binidx + i)] += 1
            if direction == "R":
                # forward direction
                for i in range(10):
                    if 50*(binidx+i) >= lengths[seq_id]:
                        continue
                    rev_scores[(seq_id, binidx + i)] += score
                    rev_counts[(seq_id, binidx + i)] += 1
            if len(fwd_scores) > 10000:
                for seq_id, binidx in fwd_scores:
                    score, count = fwd_scores[(seq_id, binidx)], fwd_counts[(seq_id, binidx)] 
                    score = score/count
                    if score > 0:
                        print(seq_id, binidx*50, binidx*50+50, score,sep="\t",file=ffwd)
                fwd_scores, fwd_counts = defaultdict(float), defaultdict(int)
            if len(rev_scores) > 10000:
                for seq_id, binidx in rev_scores:
                    score, count = rev_scores[(seq_id, binidx)], rev_counts[(seq_id, binidx)]
                    score = score/count
                    if score > 0:
                        print(seq_id, binidx*50, binidx*50+50, score, sep="\t",file=frev)
                rev_scores, rev_counts = defaultdict(float), defaultdict(int)

    for seq_id, binidx in fwd_scores:
        score, count = fwd_scores[(seq_id, binidx)], fwd_counts[(seq_id, binidx)]
        score = score/count
        if score > 0:
            print(seq_id, binidx*50, binidx*50+50, round(score,2),sep="\t",file=ffwd)

    for seq_id, binidx in rev_scores:
        score, count = rev_scores[(seq_id, binidx)], rev_counts[(seq_id, binidx)]
        score = score/count
        if score > 0:
            print(seq_id, binidx*50, binidx*50+50, round(score,2),sep="\t",file=frev)    
    ffwd.close()
    frev.close()

if __name__ == "__main__":
    main()
