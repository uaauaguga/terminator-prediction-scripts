#!/usr/bin/env python
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(description='make fake coordiante for transterm')
    parser.add_argument('--input','-i',type=str, required=True, help='input fasta file')
    parser.add_argument('--output','-o',type=str, required=True, help="output fasta file")
    args = parser.parse_args()
    fin = open(args.input)
    lengths = {} 
    for line in fin:
        if line.startswith(">"):
            attrs = line.strip()[1:].split(" ")
            seq_id = attrs[0]
            if seq_id in lengths:  #"non unique sequence id detected"
                continue
            lengths[seq_id] = 0
        else:
            lengths[seq_id] += len(line.strip())
    fout = open(args.output,"w")
    for seq_id in lengths:
        L = lengths[seq_id]
        name = f"{seq_id}-fake-1"
        fout.write(f"{name}\t1\t2\t{seq_id}\n")
        name = f"{seq_id}-fake-2"
        fout.write(f"{name}\t{L-1}\t{L}\t{seq_id}\n")
    fin.close()
    fout.close()
            

if __name__ == "__main__":
    main() 
