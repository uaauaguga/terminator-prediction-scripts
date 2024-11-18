#!/usr/bin/env python
import pandas as pd
from copy import deepcopy
from sklearn.metrics import roc_auc_score
import pickle
import numpy as np
from pyopls import OPLS, OPLSValidator
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import train_test_split

def main():
    positive = pd.read_csv("tools/mvRDT/positive.features.txt",sep="\t",index_col=0)
    negative = pd.read_csv("tools/mvRDT/background/features.txt",sep="\t",index_col=0)
    features = open("tools/mvRDT/features.txt").read().strip().split("\n")
    selected_negative = negative.iloc[:200,:]
    instances = pd.concat([positive.loc[:,features],selected_negative.loc[:,features]])
    y = np.concatenate([np.ones(positive.shape[0]),np.zeros(selected_negative.shape[0])])
    X = instances.values
    # sanity check: a rough test of fitting 
    for i in range(50):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=i)
        opls = OPLS(37)
        lda = LDA(solver='lsqr',shrinkage='auto')
        opls = opls.fit(X_train, y_train)
        Z = opls.transform(X_train)
        clf = lda.fit(Z, y_train)
        Z = opls.transform(X_test)
        y_pred = clf.predict_proba(Z)[:,1]
        AUROC = roc_auc_score(y_test, y_pred)
        print(AUROC)
    models = []
    for i in range(10):
        selected_negative = negative.iloc[i*200:(i+1)*200,:]
        instances = pd.concat([positive.loc[:,features],selected_negative.loc[:,features]])
        y = np.concatenate([np.ones(positive.shape[0]),np.zeros(selected_negative.shape[0])])
        X = instances.values
        lda = LDA(solver='lsqr',shrinkage='auto')
        opls = OPLS(37)
        opls = opls.fit(X, y)
        Z = opls.transform(X)
        clf = lda.fit(Z, y)
        model = {"opls":opls,"clf":clf}
        models.append(deepcopy(model))
    with open('tools/mvRDT/model.PLSDA.shuffled.bg.10.ensemble.pkl', 'wb') as f:
        pickle.dump(models,f)


if __name__ == "__main__":
    main()
