#!/usr/bin/env python

def main():
    # reg 128: only use 101-114 for validation, not use for training
    # or only use 0-100 for training, not use for validation
    # RF00174-Cobalamin::OTU-1427:3300027967.a:Ga0209272_10001117:4719-4901(+)	246	296	Fusobacteriota	64 
        
    fout = open("output/back-align-intervals/train.bed","w")
    with open("output/back-align-intervals/0-100.bed") as f:
        for line in f:
            # OTU-9255:NC_010002_3019:0027::OTU-9255:NC_010002:3312081-3312225(+)
            i = int(line.split("::")[0].split(":")[-1])
            if  i <= 100:
                fout.write(line)

    fout = open("output/back-align-intervals/validation.bed","w")
    with open("output/back-align-intervals/101-114.bed") as f:
        for line in f:
            # OTU-9255:NC_010002_3019:0027::OTU-9255:NC_010002:3312081-3312225(+)
            i = int(line.split("::")[0].split(":")[-1])
            if i > 100:
                fout.write(line)

if __name__ == "__main__":
    main()
