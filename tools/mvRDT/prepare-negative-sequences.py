#!/usr/bin/env python
from pyfaidx import Fasta
import numpy as np
from ushuffle import shuffle
np.random.seed(666)

def main():

    fout = open("tools/mvRDT/background/chunked.fa","w")
    n = 0
    for asm_id in ["GCF_000006945.1","GCF_000005845.2"]:
        fasta = Fasta(f"dataset/genome/refseq/assemblies-short/{asm_id}.fna")
        sequence = fasta[list(fasta.keys())[0]]
        for i in range(1000):
            p = np.random.randint(len(sequence)-500)
            s = sequence[p:p+500]
            if np.random.rand() > 0.5:
                s = s.reverse.complement
            s = str(s)
            s = shuffle(s.encode(),2).decode()
            print(">" + str(n).zfill(4), file=fout)
            print(s, file=fout)
            n += 1        
    fout.close()
            
         


if __name__ == "__main__":
    main()
