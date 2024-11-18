#!/usr/bin/env python
import argparse
import logging
import numpy as np
from pyfaidx import Fasta
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract sequences')

def main():
    parser = argparse.ArgumentParser(description='extract terminator flanking sequences')
    parser.add_argument('--bed','-b',type=str, required=True, help="annotated terminators in bed format")
    parser.add_argument('--genome','-g',type=str, required=True, help="input genome sequences")
    parser.add_argument('--fasta', '-f', type=str, required=True, help="output terminator sequences")
    parser.add_argument('--left', '-l', type=int, default=100, help="left flanking length")
    parser.add_argument('--right', '-r', type=int, default=300, help="right flanking length")
    args = parser.parse_args()

    logger.info("Load genomes ...")
    fasta = Fasta(args.genome)

    fout = open(args.fasta,"w")
    logger.info("Extract terminator sequences ...")
    #NODE_143_length_42442_cov_6.49621       8527    8585    .       0.902   -       <|<     29,8    concordant      downstream      1_7|1_8
    with open(args.bed) as f:
        for line in f:
            fields = line[:-1].split("\t")
            if fields[8] != "concordant":
                continue
            if fields[9] not in ["gene:3'","downstream"]:
                continue
            seq_id, start, end = fields[:3]
            start, end = int(start), int(end)
            strand = fields[5]
            if fields[9] == "downstream":
                if strand == "+":
                    gene_id = fields[10].split("|")[0]
                else:
                    gene_id = fields[10].split("|")[1]
            else:
                gene_id = fields[10]
            if strand == "+":
                start -= args.left
                end += args.right
            else:
                start -= args.right
                end += args.left
            if (start < 0) or (end > len(fasta[seq_id])):
                continue
            sequence = fasta[seq_id][start:end]
            if strand == "-":
                sequence = sequence.reverse.complement            
            sequence = str(sequence)            
            print(f">{seq_id}:{start}-{end}({strand}):{gene_id} {fields[4]}",file=fout) 
            print(sequence, file=fout)
    fout.close()                         

if __name__ == "__main__":
    main()
