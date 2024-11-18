#!/usr/bin/env python
from statsmodels.stats.multitest import fdrcorrection
import argparse

def load_records(path,strand):
    records = []
    with open(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, start, end, odd, pvalue = fields[:5]
            start, end, odd, pvalue = int(start), int(end), float(odd), float(pvalue)
            records.append([seq_id, start, end, odd, pvalue, strand])
    return records

def main():
    parser = argparse.ArgumentParser(description='combine regions and calculate FDR')
    parser.add_argument('--forward','-f',type=str,required=True, help='forward signal')
    parser.add_argument('--reverse','-r',type=str,required=True, help='reverse signal')
    parser.add_argument('--output','-o',type=str,
                        required=True, help = 'combined signal in bed format')
    args = parser.parse_args()
    #NC_000913.3     175     225     0.0     1.0     0       1       3       16
    
    records = load_records(args.forward,"+") + load_records(args.reverse,"-")
    
    pvalues = [r[4] for r in records]
    _, fdrs = fdrcorrection(pvalues,method='indep')
    for i in range(len(fdrs)):
        records[i].append(fdrs[i])

    fout = open(args.output,"w")
    records = sorted(records,key=lambda x:(x[0],x[1],x[2]))
    for r in records:
        print(*r,sep="\t",file=fout)
    fout.close()

        
if __name__ == "__main__":
    main()        
     

if __name__ == "__main__":
    main()
