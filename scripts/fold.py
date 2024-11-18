#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
import io
from multiprocessing import Pool
import numpy as np
from ushuffle import shuffle, Shuffler
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('calulate energy')


def energy(sequence):
    cmd = ["RNAfold"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    lines = proc.communicate(sequence.encode())[0].decode()
    line = lines.strip().split("\n")[-1]
    e = line[line.rfind("(")+1:line.rfind(")")]
    proc.wait()
    return float(e)

def main():
    parser = argparse.ArgumentParser(description='calculate folding energy')
    parser.add_argument('--input','-i',help="Input fasta")
    parser.add_argument('--output','-o',required=True,help="Output path for predicted structure in dot bracket notation")
    parser.add_argument('--jobs','-j', type=int, default=4, help = "Threads for processing.")
    global args
    args = parser.parse_args()

    sequences = {}
    logger.info("load fasta ...")
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split(" ")[0]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip()

    logger.info("calculate energy ...")
    workers = {}
    pool = Pool(args.jobs)
    seq_ids = list(sequences.keys())
    for seq_id in sequences:
        sequence = sequences[seq_id]
        workers[seq_id] = pool.apply_async(func=energy, args=(sequence,))

    fout = open(args.output,"w")
    for seq_id in workers:
        e = workers[seq_id].get()
        print(seq_id, e, len(sequences[seq_id]), sep="\t",file=fout)
        fout.flush()
    fout.close()


if __name__ == "__main__":
    main()
