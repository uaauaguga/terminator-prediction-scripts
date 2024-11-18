#!/usr/bin/env python
import argparse
import logging
import re
import os
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract hits')

def extract_hits(tbl):
    #Rho N-terminal: 24.9; Rho RNA-binding: 27; ATPase: 29.5
    cutoffs = {"Rho_N":24.9,"Rho_ATPase":29.5,"Rho_RNA_bind":27}
    domains_by_protein = {}
    with open(tbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            #3300027510.a:Ga0209537_1016864_1 -          Rho_N                PF07498.16   4.4e-10   36.3    
            fields = re.split(r"\s+",line)
            pfam_domain = fields[2]
            if pfam_domain not in cutoffs:
                continue
            score = float(fields[5])
            if score < cutoffs[pfam_domain]:
                continue
            protein_id = fields[0]            
            if protein_id not in domains_by_protein:
                domains_by_protein[protein_id] = []
            domains_by_protein[protein_id].append(pfam_domain)
    records = []
    for protein_id in domains_by_protein:
        domains = domains_by_protein[protein_id]
        domains = ",".join(sorted(list(domains)))
        records.append((protein_id,domains))
    if len(records) == 0:
        records = [("","")]
    return records

def main():
    parser = argparse.ArgumentParser(description='extract hits from hmmsearch table')
    parser.add_argument('--input-directory','-id',type=str, required=True, help="input hmmsearch table")
    parser.add_argument('--output', '-o', type=str, required=True, help="output hits profile")
    args = parser.parse_args()

    fout = open(args.output,"w")
    for tbl in sorted(os.listdir(args.input_directory),key=lambda x:x.split("-")[-1][:-4]):
        path = os.path.join(args.input_directory, tbl) 
        records = extract_hits(path)
        otu_id = tbl[:-4]
        for protein_id, domains in records:
            print(otu_id, protein_id, domains, sep="\t", file=fout)
    fout.close()                 

if __name__ == "__main__":
    main()
