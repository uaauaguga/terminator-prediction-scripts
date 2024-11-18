#!/usr/bin/env python
import argparse
import os
import sys
import re
import logging
from multiprocessing import Pool
import subprocess
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('search with ensemble of cm models')



def search(motif_id):
    tbl = os.path.join(args.output_directory , motif_id + ".tbl")
    msa = os.path.join(args.output_directory , motif_id + ".stk")
    fout = open(os.path.join(args.output_directory , motif_id + ".txt"),"w")
    cm = os.path.join(args.motif_directory, motif_id + ".cm")
    logger.info(f"run cmsearch for {motif_id}...")
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/cmsearch", "-A", msa,  "--tblout", tbl,
          "--noali", "--cpu", "6", "--toponly", cm, args.fasta ]    
    p = subprocess.run(cmd,stderr=subprocess.STDOUT,stdout=fout)
    logger.info(f"convert {motif_id} hit in tabular format to bed format ...")
    bed = os.path.join(args.output_directory , motif_id + ".bed")
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/python",
           "scripts/rfam-tbl2bed.py","-i",tbl,"-o",bed]
    p = subprocess.run(cmd,stderr=subprocess.STDOUT,stdout=fout)
    fout.close()
    logger.info(f"finish processing {motif_id} .")
    return p.returncode 



def main():
    global args
    parser = argparse.ArgumentParser(description='run cmfinder pipeline')
    parser.add_argument('--motif-directory','-m',type=str, required=True, help='directory contains motif cm models')
    parser.add_argument('--fasta','-f',type=str, required=True, help='input fasta')
    parser.add_argument('--output-directory','-o',type=str, required=True, help='output directory')
    parser.add_argument('--jobs','-j',type=int,default=16,help="number of jobs to run")
    args = parser.parse_args()

    logger.info(f"use motifs in {args.motif_directory} to search against {args.fasta} ...")
    logger.info(f"run {args.jobs} jobs for cmsearch .")
    
    pool = Pool(args.jobs)
    workers = []
    for motif in os.listdir(args.motif_directory):
        if not motif.endswith(".cm"):
            continue
        motif_id = motif[:-3]
        workers.append(pool.apply_async(func=search, args=(motif_id,)))
    for worker in workers:
        code = worker.get()
    logger.info("All done .")
    
    
                 
    
    

if __name__ == "__main__":
    main()
