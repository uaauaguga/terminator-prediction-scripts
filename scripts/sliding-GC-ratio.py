#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('G/C slider')
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="calculate G/C ratio in a sliding window")
    parser.add_argument('--input','-i', help="Input fasta file", required=True)
    parser.add_argument('--output','-o', help="Output bed file",required=True)
    parser.add_argument('--window', '-w',default=72,type=int)
    parser.add_argument('--stride', '-s',default=10,type=int)
    args = parser.parse_args()
    sequences = {}
    logger.info("load sequences ...")
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split(" ")[0]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip()
    fout = open(args.output,"w")
    for seq_id in sequences:
        sequence = sequences[seq_id]
        logger.info(f"calculate GC ratio for {seq_id} ...")
        p = int(args.window/2) 
        last_p = len(sequence) - int(args.window/2)
        while p < last_p:
            start, end = p - int(args.window/2), p + int(args.window/2)
            s, e = p - int(args.stride/2), p + int(args.stride/2)
            segment = sequence[start:end]
            A, C, G, T = segment.count("A"),segment.count("C"), segment.count("G"), segment.count("T")
            A, C, G, T = 0.5 + float(A), 0.5 + float(C), 0.5 + float(G), 0.5 + float(T)
            total = A + C + G + T
            score = np.round(np.log(C/G),3)
            G_ratio = np.round(G/total,3)
            C_ratio = np.round(C/total,3)
            print(seq_id, s, e, score, G_ratio, C_ratio,sep="\t",file=fout)
            p += args.stride 
    fout.close()


if __name__ == "__main__":
    main()
