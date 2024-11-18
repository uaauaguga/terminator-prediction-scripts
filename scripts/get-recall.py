#!/usr/bin/env python
import argparse
import logging
import numpy as np
import io
import subprocess
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("calculate recall")

def get_scores(true_path, pred_path, stranded = False):
    logger.info(f"annotate {true_path} with prediction {pred_path} ...")
    cmd = ["bedtools","intersect","-loj","-a",true_path,"-b",pred_path]
    if stranded:
        cmd += ["-s"]
    print(" ".join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    scores = []
    for line in io.TextIOWrapper(p.stdout, encoding="utf-8"):
        fields = line.strip().split("\t")
        try:
            score = float(fields[10])
        except:
            score = 0
        scores.append(score)
    code = p.poll()
    return scores


def main():
    parser = argparse.ArgumentParser(description='calculate recall given positive instance and predictions')
    parser.add_argument('--positive', '-p', type=str, required=True, help='ground truth intervals in bed format')
    parser.add_argument('--terminator', '-t', type=str, required=True, help='intervals of predicted terminators in bed format, a score field should present')
    parser.add_argument('--output', '-o', type=str, required=True, help='where to save the output')
    parser.add_argument('--bins', '-b', type=int, default=50, help='number of bins to use')
    parser.add_argument('--min-score', '-ms', type=float, default=50, help='set min score to this value')
    parser.add_argument('--stranded', '-s', action= "store_true", help='whether consider strandness')
    args = parser.parse_args()

    logger.info("assign predicted scores to  ground truth regions ...")
    scores = get_scores(args.positive, args.terminator, args.stranded)
    n = len(scores)
    scores = np.array(scores)
    scores[scores < args.min_score] = args.min_score
    counts, scores = np.histogram(scores,bins=args.bins)
    counts = counts[::-1].cumsum()[::-1] 
    logger.info("saving results ...")
    with open(args.output,"w") as fout:
        for c, s in zip(counts, scores[:-1]):
            print(s,c/n,file=fout,sep="\t")
    logger.info("all done .")

if __name__ == "__main__":
    main()
