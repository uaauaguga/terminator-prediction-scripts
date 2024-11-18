#!/usr/bin/env python
import pandas as pd
def main():
    table = pd.read_csv("tools/mvRDT/2018.NAR.E.coli.train.txt",sep="\t",index_col="LabName")
    labels = table["Signal"].values
    positive_ids = table[labels != "Nothing"].index
    negative_ids = table[labels == "Nothing"].index
    
    table2 = pd.read_csv("tools/mvRDT/train/features.txt",sep="\t",index_col="names-seq") 

    features = open("tools/mvRDT/features.txt").read().strip().split("\n")
    positive_features = table2[table2.index.isin(positive_ids)].loc[:,features]
    negative_features = table2[table2.index.isin(negative_ids)].loc[:,features]

    positive_features.to_csv("tools/mvRDT/positive.features.txt",sep="\t")
    negative_features.to_csv("tools/mvRDT/validated.negative.features.txt",sep="\t")

if __name__ == "__main__":
    main()
