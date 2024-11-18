#!/usr/bin/env python
import subprocess 
import argparse
from pyfaidx import Fasta
import subprocess

def energy(sequence):
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/RNAfold"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    lines = proc.communicate(sequence.encode())[0].decode()
    line = lines.strip().split("\n")[-1]
    energy = line[line.rfind("(")+1:line.rfind(")")]
    proc.wait()
    try:
        energy = float(energy)
    except:
        print(energy)
        energy = 0
    return energy

def main():
    parser = argparse.ArgumentParser(description='evaluate folding energy near termination sites')
    parser.add_argument('--input','-i',type=str,required=True, help="3' ends in bed format")
    parser.add_argument('--fasta','-f',type=str,required=True, help='genome sequence')
    parser.add_argument('--output','-o',type=str,required=True, help = 'output sites with folding energy annotation')
    args = parser.parse_args()

    fasta = Fasta(args.fasta)
    i = 0
    fout = open(args.output,"w")
    with open(args.input) as f:
        for line in f:
            if i%100 == 0:
                print(f"{int(i/100)/10} sites processed .")
            i += 1
            seq_id, start, end, name, score, strand = line.strip().split("\t")[:6]
            start, end = int(start), int(end)
            if strand == "+":
                start, end = end - 45, end + 5
            else:
                start, end = start - 5, start + 45
                sequence = fasta[seq_id][max(start-20,0):min(start+20,len(fasta[seq_id]))]
            if start < 0 or end > len(fasta[seq_id]):
                continue
            sequence = fasta[seq_id][start:end]
            if strand == "-":
                sequence = sequence.reverse.complement
            sequence = str(sequence).replace("T","U")
            mfe = energy(sequence)
            print(line.strip(),mfe,sep="\t",file=fout)
    fout.close()
            
            

        
if __name__ == "__main__":
    main()        
     

if __name__ == "__main__":
    main()
