#!/usr/bin/env python
import os
from collections import defaultdict
import argparse 
import re
import numpy as np
from tqdm import tqdm
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("get picked motifs")
from Bio import AlignIO 
from shutil import copy

def main():
    parser = argparse.ArgumentParser(description='get selected motifs')
    parser.add_argument('--stockholm-directory','-sd',type=str, required=True, help='directory contains stockholm file')
    parser.add_argument('--cm-directory','-cd',type=str, required=True, help='input directory contains covariance model')
    parser.add_argument('--output-directory','-o',type=str, required=True, help='output directory')
    parser.add_argument('--hit-gene-number','-n',type=str, required=True, help='hitted genes')
    args = parser.parse_args()


    counts = {}
    indices = {}
    with open(args.hit_gene_number) as f:
        # class-ABY1	1	1	106
        for line in f:
            clade_id, loop, index, count = line.strip().split("\t") 
            count = int(count)
            if (clade_id, loop) not in counts:
                counts[(clade_id, loop)] = count
                indices[(clade_id, loop)] = index
            elif count > counts[(clade_id, loop)]:
                counts[(clade_id, loop)] = count
                indices[(clade_id, loop)] = index

    for clade_id, loop in tqdm(counts):
        count = counts[(clade_id, loop)]
        index = indices[(clade_id, loop)]
        if count >= 16: 
            stockholm_path = f"{args.stockholm_directory}/{clade_id}.fa.motif.h{loop}_{index}"
            if not os.path.exists(f"{args.output_directory}/{clade_id}-{loop}-{index}.stk"):
                copy(stockholm_path,f"{args.output_directory}/{clade_id}-{loop}-{index}.stk")
            cm_path = f"{args.cm_directory}/{clade_id}-{loop}-{index}.cm"
            if not os.path.exists(f"{args.output_directory}/{clade_id}-{loop}-{index}.cm"):
                copy(cm_path,f"{args.output_directory}/{clade_id}-{loop}-{index}.cm")


if __name__ == "__main__":
    main()
