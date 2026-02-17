#!/usr/bin/env python
import argparse
import numpy as np
from tqdm import tqdm
import HTSeq
import logging

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("5' end coverage")


def main():
    parser = argparse.ArgumentParser(description="get coverage of reads 5' ends")
    parser.add_argument('--input','-i',required=True,help="Input bam file at genome coordiante")
    parser.add_argument('--forward','-f',required=True,help="Output coverage in forward strand, bedgraph format")
    parser.add_argument('--reverse','-r',required=True,help="Output coverage in reverse strand, bedgraph format")
    args = parser.parse_args() 
    cvg = HTSeq.GenomicArray("auto", stranded=True, typecode='i')            
    bam_reader = HTSeq.BAM_Reader(args.input)
    n_total = 0
    logger.info("Load 5' coverage ...")
    for aln in tqdm(bam_reader):
        if aln.iv.strand == "+":
            five_prime_end = HTSeq.GenomicInterval(aln.iv.chrom, aln.iv.start , aln.iv.start + 1, aln.iv.strand)
        else:
            five_prime_end = HTSeq.GenomicInterval(aln.iv.chrom, aln.iv.end - 1, aln.iv.end, aln.iv.strand)
        cvg[five_prime_end] += 1
        n_total += 1
    logger.info(f"{n_total} reads processed .")
    logger.info("Saving coverages in forward strand ...")
    cvg.write_bedgraph_file(args.forward, strand="+")
    logger.info("Saving coverages in reverse strand ...")
    cvg.write_bedgraph_file(args.reverse, strand="-")
    logger.info("All done .")
                       
        
if __name__ == "__main__":
    main()

