#!/usr/bin/env python
import argparse
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("pick window")

def pick_local_max(scores):
    scores = np.array(scores)
    L = scores.shape[0]
    shift_left = np.concatenate([scores,[0]])[1:]
    shift_right = np.concatenate([[0],scores])[:-1]
    mask = (scores > shift_left) & (scores > shift_right)
    local_max_indices = np.where(mask)[0] 
    return local_max_indices


def main():
    parser = argparse.ArgumentParser(description='pick local max')
    parser.add_argument('--forward', '-f', type=str, required=True, help='forward scores')
    parser.add_argument('--reverse', '-r', type=str, required=True, help='reverse scores')
    parser.add_argument('--output', '-o', type=str, required=True, help='output bed file')
    args = parser.parse_args()

    logger.info("load forward strand data ...")
    fwd_locations = {}
    fwd_scores = {}
    with open(args.forward) as f:
        for line in f:
            seq_id, start, end, score = line.strip().split("\t")
            start, end, score = int(start), int(end), float(score)
            if seq_id not in fwd_locations:
                fwd_locations[seq_id] = []
                fwd_scores[seq_id] = []
            fwd_locations[seq_id].append((start, end))
            fwd_scores[seq_id].append(score)

    logger.info("load reverse strand data ...")
    rev_locations = {}
    rev_scores = {}
    with open(args.reverse) as f:
        for line in f:
            seq_id, start, end, score = line.strip().split("\t")
            start, end, score = int(start), int(end), float(score)
            if seq_id not in rev_locations:
                rev_locations[seq_id] = []
                rev_scores[seq_id] = []
            rev_locations[seq_id].append((start, end))
            rev_scores[seq_id].append(score)

    seq_ids = set()
    for seq_id in fwd_locations:
        seq_ids.add(seq_id)
    for seq_id in rev_locations:
        seq_ids.add(seq_id)
    seq_ids = sorted(list(seq_ids))      
    logger.info("pick local max ...")
    logger.info(f"results will be saved to {args.output} ...")
    fout = open(args.output,"w")
    for seq_id in seq_ids:
        print(seq_id)
        records = []
        if seq_id in fwd_scores:
            for i in pick_local_max(fwd_scores[seq_id]):
                start, end = fwd_locations[seq_id][i]
                score = fwd_scores[seq_id][i]
                records.append((seq_id, start, end, ".", score, "+"))
        if seq_id in rev_scores:
            for i in pick_local_max(rev_scores[seq_id]):
                start, end = rev_locations[seq_id][i]
                score = rev_scores[seq_id][i]
                records.append((seq_id, start, end, ".", score, "-"))
        records = sorted(records,key=lambda x:(x[1],x[2]))
        for fields in records:
            print(*fields,sep="\t",file=fout)
    fout.close()
         

        
if __name__ == "__main__":
    main()
