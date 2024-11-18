#!/usr/bin/env python
import argparse
import numpy as np
from scipy.stats import fisher_exact
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('rhoterm-regions')

def load_coverages(path):
    coverages = {}
    with open(path) as f:
        for line in f:
            seq_id, start, end, count = line.strip().split("\t")
            count = int(count)
            if seq_id not in coverages:
                coverages[seq_id] = []
            coverages[seq_id].append(count)
    return coverages

def merge_coverage(coverages):
    combined = {}
    for contig_id in coverages[0]:
        combined[contig_id] = np.zeros(len(coverages[0][contig_id]),dtype=int)
    for c in coverages:
        for contig_id in c:
            combined[contig_id] += np.array(c[contig_id])
    return combined


def main():
    parser = argparse.ArgumentParser(description='identify rho termination region')
    parser.add_argument('--treatments','-ts',type=str, required=True, help='coverage files of BCM treatment/rho KO, separated by comma')
    parser.add_argument('--controls','-cs',type=str, required=True, help="coverage files of controls, separated by comma")
    parser.add_argument('--output','-o',type=str, required=True, help="putative rhoterm regions")
    parser.add_argument('--minimal-counts','-mc',type=int, default=0, help="minimal read counts of two bins to consider")
    parser.add_argument('--window','-w',type=int, default=50, help="bin size")
    parser.add_argument('--bins','-b',type=int, default=4, help="number of bins to consider at each side")
    parser.add_argument('--strand','-s',type=str, choices = ["+","-"],required=True, help="considered strand")
    args = parser.parse_args()

    logger.info("load treatment data ...")
    treatment_coverages = []
    for path in args.treatments.split(","):
        logger.info(f"load coverage from {path} ...")
        treatment_coverages.append(load_coverages(path))
    treatment_coverages = merge_coverage(treatment_coverages)

    logger.info("load control data ...")
    control_coverages = []
    for path in args.controls.split(","):
        logger.info(f"load coverage from {path} ...")
        control_coverages.append(load_coverages(path))
    control_coverages = merge_coverage(control_coverages)

    fout = open(args.output,"w")
    logger.info("perform fisher exact test ...")
    contig_ids = [contig_id for contig_id in treatment_coverages if contig_id in control_coverages]
    for contig_id in contig_ids:
        nbins = len(control_coverages[contig_id])
        logger.info(f"processing {contig_id} ...")
        for i in range(args.bins,nbins-args.bins):
            if (i+1)%5000 == 0:
                logger.info(f"{int((i+1)/1000)} K bins processed .")
            t1, t2 = treatment_coverages[contig_id][i-args.bins:i].sum(), treatment_coverages[contig_id][i:i+args.bins].sum()
            c1, c2 =   control_coverages[contig_id][i-args.bins:i].sum(),   control_coverages[contig_id][i:i+args.bins].sum()
            if (t1 + t2) < args.minimal_counts or (c1 + c2) < args.minimal_counts:
                odd, pvalue = np.nan, 1
            else:
                if args.strand == "+":
                    odd, pvalue = fisher_exact([[t2, t1],[c2, c1]],alternative="greater")
                else:
                    odd, pvalue = fisher_exact([[t1, t2],[c1, c2]],alternative="greater")
            start, end = int(args.window*(i-0.5)), int(args.window*(i+0.5))
            print(contig_id, start, end, odd, pvalue, t1, t2, c1, c2, sep="\t",file=fout)
            fout.flush()

    logger.info("all done .")
        

if __name__ == "__main__":
    main()
