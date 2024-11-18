#!/usr/bin/env python
import argparse
import os
import sys
from tqdm import tqdm
import logging
from multiprocessing import Pool
from glob import glob
import subprocess
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('run cmfinder')


cmfinder_exe = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/cmfinder-env/bin/cmfinder04.pl" #"~/qhsky1/miniconda/envs/cmfinder-env/bin/cmfinder04.pl"


def cmfinder(fasta_path,log_path):
    flog = open(log_path,"wb")
    cwd = os.getcwd()
    fasta = os.path.basename(fasta_path)
    indir = os.path.dirname(fasta_path)
    os.chdir(indir)
    cmd = [cmfinder_exe,"-f","0.5",fasta]
    logger.info(f"run cmfinder for {fasta_path} ...")
    p = subprocess.run(cmd,stderr=subprocess.STDOUT,stdout=flog)
    flog.close()
    os.chdir(cwd)
    if p.returncode == 0:
        logger.info(f"cmfinder run successfully for {fasta_path} .")
        with open(fasta_path + ".checkpoint","w") as f:
            f.write("done")
    else:
        logger.info(f"cmfinder failed for {fasta_path} .")
        with open(fasta_path + ".checkpoint","w") as f:
            f.write("failed")
    return int(p.returncode)
    


def main():
    parser = argparse.ArgumentParser(description='run cmfinder pipeline')
    parser.add_argument('--input-directory','-i',type=str, required=True, help='input directory')
    parser.add_argument('--log-directory','-l',type=str,help="log directory")
    parser.add_argument('--jobs','-j',type=int,default=16,help="number of jobs to run")
    parser.add_argument('--dry-run',action="store_true",help="only determine number of jobs to run")
    args = parser.parse_args()

    logger.info("get sequence path ...") 
    fasta_paths = []
    log_paths = []

    processed_clade_ids = set()
    clade_ids = set()
    for file_name in os.listdir(args.input_directory):
        if "fa.motif.h" in file_name:
            clade_id = file_name[:-11]        
            processed_clade_ids.add(clade_id)
            continue
        if not file_name.endswith(".fa"):
            continue
        clade_id = file_name[:-3]
        clade_ids.add(clade_id)
    n_species_level = 0
    n_processed = len(processed_clade_ids)
    for clade_id in clade_ids:
        if clade_id.startswith("species"):
            n_species_level += 1
            continue
        if clade_id in processed_clade_ids:
            continue
        if os.path.exists(os.path.join(args.input_directory,clade_id + ".fa.checkpoint")):
            n_processed += 1
            continue
        fasta_paths.append(os.path.join(args.input_directory,clade_id + ".fa"))
        log_paths.append(os.path.join(args.log_directory,clade_id + ".log"))

    logger.info(f"{n_processed} clades already processed .")
    logger.info(f"skip {n_species_level} species level clades .")
    logger.info(f"{len(fasta_paths)} clades to process .")

    if args.dry_run:
        sys.exit(0)

    if len(fasta_paths) == 0:
        logger.info(f"no fasta file present in {args.input_directory}, exiting .")
        sys.exit(1) 

    
    pool = Pool(args.jobs)
    workers = []
    for fasta_path, log_path in zip(fasta_paths,log_paths):
        workers.append(pool.apply_async(func=cmfinder, args=(fasta_path,log_path,)))
    for worker in workers:
        code = worker.get()
    
    logger.info(f"All done.")                 
    
    

if __name__ == "__main__":
    main()
