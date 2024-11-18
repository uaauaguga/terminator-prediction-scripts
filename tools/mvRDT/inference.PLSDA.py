#!/usr/bin/env python
import argparse
from lightgbm import LGBMClassifier
import pickle
import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='run predictions')
    parser.add_argument('--input', '-i', type=str, required=True, help='features')
    parser.add_argument('--output', '-o', type=str, required=True, help='scores')
    parser.add_argument('--batch-size', '-bs', type=int, default=4096, help='batch size for inference')
    parser.add_argument('--model', '-m', type=str, required=True, help='model used for inference')
    args = parser.parse_args()
    with open(args.model,"rb") as f:
        models = pickle.load(f)
    table = pd.read_csv(args.input,sep="\t",index_col=0)
    features = open("tools/mvRDT/features.txt").read().strip().split("\n")
    cursor = 0
    fout = open(args.output,"w")
    while cursor < table.shape[0]:        
        X = table.iloc[cursor:cursor+args.batch_size,:].loc[:,features].values
        seq_ids = table.index[cursor:cursor+args.batch_size]        
        y_pred_voted = np.zeros(X.shape[0])
        n = 0
        for model in models:
            opls = model["opls"]
            clf = model["clf"]
            Z = opls.transform(X)
            y_pred = clf.predict_proba(Z)[:,1]
            y_pred_voted += y_pred
            n += 1
        y_pred_voted = y_pred_voted/n
        probs = np.round(y_pred_voted,2)
        for seq_id, prob in zip(seq_ids, probs):
            print(seq_id, prob, sep="\t", file=fout)
        cursor += args.batch_size
    fout.close()
    
             
if __name__ == "__main__":
    main()
