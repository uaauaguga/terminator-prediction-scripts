#!/usr/bin/env python
import argparse
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("pick best bin from overlapping ones")

from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='pick best interval from overlapped ones')
    parser.add_argument('--input', '-i',required=True,help="input intervals")
    parser.add_argument('--output','-o',required=True,help="output intervals")
    args = parser.parse_args()


    logger.info(f"load intervals from {args.input} ...")
    logger.info(f"picked intervals will be saved to {args.output} .")
    fin = open(args.input)
    fout = open(args.output,"w")    
    starts, ends = defaultdict(list), defaultdict(list)
    strands = defaultdict(list)
    scores = defaultdict(list)
    names = defaultdict(list)
    for line in fin:
        fields = line.strip().split("\t")
        seq_id, start, end, name, score, strand = fields[:6]
        start, end, score = int(start), int(end), float(score)
        starts[(seq_id,strand)].append(start)
        ends[(seq_id,strand)].append(end)
        scores[(seq_id,strand)].append(score)
        names[(seq_id,strand)].append(name)
    #lengths = {} 
    #for seq_id in ends:
    #    lengths[seq_id] = max(ends[seq_id]) 
    
    for seq_id, strand in scores:
        values = np.array(scores[(seq_id,strand)])
        ninty = sorted(values)[int(0.95*len(values))]
        shift_right = np.array([0] + scores[(seq_id,strand)])
        shift_left = np.array(scores[(seq_id,strand)] + [0])
        mask = (values > shift_right[:-1]) & (values > shift_left[1:]) & (values >ninty)
        print(mask.sum())
        indices = np.where(mask)[0]
        for i in indices:
            start, end, name, score= starts[(seq_id,strand)][i], ends[(seq_id,strand)][i], names[(seq_id,strand)][i], scores[(seq_id,strand)][i]
            print(seq_id, start, end,name,score,strand,sep="\t",file=fout)
    fin.close()
    fout.close()

if __name__ == "__main__":
    main()

