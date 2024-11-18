#!/usr/bin/env python
import argparse
from scipy.interpolate import interp1d
import numpy as np


def load_scores(path):
    cutoffs, scores = [], []
    with open(path) as f:
        for line in f:
            cutoff, score = line.strip().split("\t")
            cutoff, score = float(cutoff), float(score)
            cutoffs.append(cutoff)
            scores.append(score)
    cutoffs, scores = np.array(cutoffs), np.array(scores)
    mask = scores!=1
    return cutoffs[mask], scores[mask]



def main():
    parser = argparse.ArgumentParser(description='get ROC curve')
    parser.add_argument('--tpr', '-t', type=str, required = True, help='true negative rate with score cutoff')
    parser.add_argument('--fpr', '-f', type=str, required = True, help='true negative with score cutoff')
    parser.add_argument('--bins', '-b', type=int, default=50, help='number of bins to use')
    parser.add_argument('--output', '-o', type=str, required = True, help='roc curve')
    args = parser.parse_args()

    fpr_cutoffs, fprs = load_scores(args.fpr)
    tpr_cutoffs, tprs = load_scores(args.tpr)

    fprinterp = interp1d(fpr_cutoffs, fprs)
    tprinterp = interp1d(tpr_cutoffs, tprs)
    min_score = max(fpr_cutoffs.min(),tpr_cutoffs.min())
    max_score = min(fpr_cutoffs.max(),tpr_cutoffs.max())

    fout = open(args.output,"w")
    for cutoff in np.linspace(min_score,max_score,args.bins):
        try:
            tpr = tprinterp(cutoff)
            fpr = fprinterp(cutoff)
            print(fpr,tpr,cutoff,sep="\t",file=fout)
        except:
            continue
    fout.close()
if __name__ == "__main__":
    main() 
