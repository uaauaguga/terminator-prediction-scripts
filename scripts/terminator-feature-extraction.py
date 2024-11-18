#!/usr/bin/env python
from tqdm import tqdm
import argparse
import re

def get_stem_loop_size(structure):
    intervals = []
    kept_stem_size = 0
    kept_loop_size = 0
    kept_interval = (-1,-1)
    for m in re.finditer(r"\(\.+\)",structure):
        l, r = m.start(),m.end()
        loop_size = r - l - 2
        left_offset = 0
        right_offset = 0
        left_consecutive_unpair, right_consecutive_unpair = 0, 0
        while True:
            while (l-left_offset>=0) and (structure[l-left_offset] == "."):
                left_consecutive_unpair += 1
                left_offset += 1
            while (r+right_offset <= len(structure)) and (structure[r+right_offset-1] == "."):
                right_consecutive_unpair += 1
                right_offset += 1
            if (left_consecutive_unpair >= 4) or (right_consecutive_unpair >= 4):
                left_offset -= left_consecutive_unpair 
                right_offset -= right_consecutive_unpair 
                break            
            if (l-left_offset< 0) or (r+right_offset > len(structure)):
                break
            if not ((structure[l-left_offset] == "(") and (structure[r+right_offset-1] == ")")):
                break
            else:
                left_consecutive_unpair = right_consecutive_unpair = 0
            left_offset += 1
            right_offset += 1
            if not ((l-left_offset>=0) and (r+right_offset-1 < len(structure))):
                break 
        substructure = structure[l-left_offset:r+right_offset]
        left_offset -= substructure.find("(")
        right_offset -= len(substructure) - substructure.rfind(")") - 1       
        stem_size =  1 + int((left_offset+right_offset)/2)
        if stem_size > kept_stem_size:
            kept_stem_size = stem_size
            kept_loop_size = loop_size
            kept_interval = (l-left_offset,r+right_offset)
    return kept_stem_size, kept_loop_size

def get_U_run_length(sequence):
    L = 0
    pattern = r"(T|t|U|u)+[^(T|t|U|u)]?(T|t|U|u)+[^(T|t|U|u)]?(T|t|U|u)+"
    #pattern = r"(T|t|U|u)+"
    for m in re.finditer(pattern,sequence):
        l = m.end() - m.start() 
        #print(sequence[m.start():m.end()])
        if l > L:
            L = l
    return L


def main():
    parser = argparse.ArgumentParser(description="refine candidate terminator")
    parser.add_argument('--input', '-i', help="candidate terminator sequence with predicted structure", required=True)
    parser.add_argument('--output','-o', help="filter and refinement",required=True)
    args = parser.parse_args()
    
    fout = open(args.output,"w") 
    print("seq id","stem length","loop length","U count", sep="\t",file=fout)
    with open(args.input) as f:
        for line in f:
            seq_id = line.strip()[1:].split(" ")[0]
            sequence = next(f).strip()
            structure = next(f).strip()
            fields = structure.split(" ")
            structure, energy = fields[0], fields[-1]
            energy = energy.replace("(","").replace(")","").strip()
            U_count = get_U_run_length(sequence)
            stem_size, loop_size = get_stem_loop_size(structure)
            print(seq_id, stem_size, loop_size, U_count, sep="\t",file=fout)
    fout.close()

if __name__ == "__main__":
    main()
