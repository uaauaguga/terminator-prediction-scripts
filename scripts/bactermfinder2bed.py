#!/usr/bin/env python
import argparse
import logging
import re

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('bactermfinder')
def main():
    parser = argparse.ArgumentParser(description='reformat')
    parser.add_argument('--input','-i',type=str,required=True, help="input prediction")
    parser.add_argument('--output','-o',type=str,required=True, help = 'output in bed format')
    args = parser.parse_args()

    fin = open(args.input)
    fout = open(args.output,"w")
    pattern = r'(\d+)_(.+)_(\d+)_(\d+)_([+-])'
    _ = next(fin)
    for line in fin:
        fields = line.strip().split(",")
        #1_NC_000913.3_0_101_+
        score = fields[5]
        sample_id = fields[0]
        #('NC_000913.3', '4641537', '4641638', '+')
        init_id, seq_id, start, end, strand = re.match(pattern,sample_id).groups()
        init_id = init_id.zfill(12)
        print(seq_id, start, end, init_id, score, strand, file=fout,sep="\t")
    fout.close()
    fin.close()

if __name__ == "__main__":
    main()
