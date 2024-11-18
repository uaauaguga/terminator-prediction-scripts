#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('chunkify sequence')

def main():
    parser = argparse.ArgumentParser(description='chunkify sequence for mvRDT')
    parser.add_argument('--input', '-i', type=str, required=True, help='input genome sequence')
    parser.add_argument('--output', '-o', type=str, required=True, help='output chunked sequences')
    args = parser.parse_args()
    rc_lut = {"A":"T","C":"G","G":"C","T":"A"}


    logger.info("load sequences ...")
    sequences = {}
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line.strip()[1:].split(" ")[0]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip()

    logger.info("save chunkified sequences ...")
    fout = open(args.output,"w")
    for seq_id in sequences:
        p = 0
        n = 1
        while p < len(sequences[seq_id])-500:
            s, e = p, p + 500
            sequence = sequences[seq_id][s:e]
            sequence_rc = "".join([ rc_lut.get(c,"N") for c in sequence ][::-1])
            fout.write(f">Seq{n}_{seq_id}_{s+1}_{e}_F\n")
            fout.write(sequence + "\n")
            n += 1
            fout.write(f">Seq{n}_{seq_id}_{s+1}_{e}_R\n")
            fout.write(sequence_rc + "\n")
            n += 1
            p += 50 
    fout.close()


if __name__ == "__main__":
    main()
