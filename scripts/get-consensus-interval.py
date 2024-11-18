#!/usr/bin/env python

import argparse
import numpy as np
from tqdm import tqdm

def rle(inarray):
    """ run length encoding. Partial credit to R rle function. 
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values) """
    ia = np.asarray(inarray)                # force numpy
    n = len(ia)
    if n == 0: 
        return (None, None, None)
    else:
        y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element posi
        z = np.diff(np.append(-1, i))       # run lengths
        p = np.cumsum(np.append(0, z))[:-1] # positions
    return(z, p, ia[i])


def main():
    parser = argparse.ArgumentParser(description="get consensus interval from hit coverage")
    parser.add_argument('--input','-i', help="input coverage", required=True)
    parser.add_argument('--output','-o', help="output in bed format",required=True)
    parser.add_argument('--count-cutoff', '-cc', type=int, default=5)
    parser.add_argument('--fraction-cutoff', '-fc', type=float, default=0.25)
    parser.add_argument('--length-cutoff', '-lc', type=int, default=16)
    args = parser.parse_args()
    fout = open(args.output,"w")
    with open(args.input) as f:
        for line in tqdm(f):
            seq_id, count, coverage = line.strip().split("\t")[:3]
            count = int(count)
            if count < args.count_cutoff:
                continue
            coverage = np.array(coverage.split(",")).astype(int)
            mask = coverage >= args.fraction_cutoff*count
            lengths, positions, values = rle(mask)      
            for position, length, value in zip(positions,lengths,values):
                if value and length > args.length_cutoff:
                    local_max = coverage[position:position+length].max()
                    fout.write(f"{seq_id}\t{position}\t{position+length}\t{local_max}\n")
    fout.close()
if __name__ == "__main__":
    main()
