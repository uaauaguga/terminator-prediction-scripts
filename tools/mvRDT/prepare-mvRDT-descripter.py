#!/usr/bin/env python
import numpy as np
import argparse
import os
import tempfile
import subprocess
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("extract descriptor")

def main():
    parser = argparse.ArgumentParser(description='generate descripter for prediction')
    parser.add_argument('--input','-i',type=str, required=True, help='input fasta file')
    parser.add_argument('--output-directory','-od',type=str, required=True, help="output bed file")
    parser.add_argument('--chunked','-c', action = "store_true", help="whether the input sequences are 500 nt chunks")
    args = parser.parse_args()
    
    if not os.path.exists(args.output_directory):
        logger.info(f"{args.output_directory} does not exists, create it")
        os.mkdir(args.output_directory)    

    if args.chunked:
        fasta = args.input
        logger.info("input sequence are already chunkified .")
    else:
        logger.info("chunkify input sequences ...")
        fasta = os.path.join(args.output_directory,"chunked.fa")
        cmd = ["scripts/prepare-chunks-for-mvRDT.py","-i", args.input, "-o", fasta]
        print(" ".join(cmd))
        subprocess.run(cmd)
    
    tns = tempfile._get_candidate_names()
    logger.info("make descriptor ...")
    tmpname = next(tns)
    cmd = ["python", "tools/mvRDT/MakeDescriptors.py", fasta , tmpname ]
    subprocess.run(cmd)    
    name = fasta.split("/")[-1][:-6]
    path = f"{name}_{tmpname}_VarOPLSDA.txt" 
    tgt = os.path.join(args.output_directory,"VarOPLSDA.txt")
    os.rename(path, tgt)
    
    logger.info("run RNA folding ...")
    energy = os.path.join(args.output_directory,"energy.txt")
    cmd = ["scripts/fold.py", "-i", fasta, "-o", energy ]     
    subprocess.run(cmd)

    logger.info("combine results ...")
    table = pd.read_csv(tgt, sep="\t",dtype={'names-seq': str})
    table = table.set_index('names-seq')
    mfe_by_kb = {}
    with open(energy) as f:
        for line in f:
            seq_id, energy, length = line.strip().split("\t")
            energy, length = float(energy), int(length)
            mfe_by_kb[seq_id] = energy*1000/length

    table["TempDenDG"] = table.index.map(lambda x:mfe_by_kb[x])
    table.to_csv(os.path.join(args.output_directory,"features.txt"),sep="\t")
     
            
    
    
    
 

if __name__ == "__main__":
    main()

