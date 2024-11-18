#!/usr/bin/env python
import argparse
import numpy as np
import pysam
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('get coverage by bins')

def write_bedgraph(coverages, path, scaler = None):
    fout = open(path,"w")
    for contig_id in sorted(list(coverages.keys())):
        nbins = coverages[contig_id].shape[0]
        last_count = 0
        start, end = 0, 0
        for i in range(nbins):            
            count = coverages[contig_id][i]
            if count == 0:
                last_count = count
                continue
            if count == last_count:
                end = i + 1
            else:
                if last_count > 0:                    
                    if scaler != None:
                        value = last_count*scaler
                    else:
                        value = last_count
                    print(contig_id, start, end, value, file=fout, sep="\t")
                start, end = i, i + 1
                last_count = count
    fout.close()


def main():
    parser = argparse.ArgumentParser(description='get coverage by bins')
    parser.add_argument('--input','-i',type=str, required=True, help='input bam file')
    parser.add_argument('--forward','-f',type=str, required=True, help="output forward coverage file")
    parser.add_argument('--reverse','-r',type=str, required=True, help="output reverse coverage file")
    parser.add_argument('--strand','-s',type=str, default="forward", choices=["forward","reverse"],help="libtype of the bam file")
    parser.add_argument('--normalize','-n',action="store_true",help="whether normalize the coverage")
    args = parser.parse_args()

    forward_coverages = {}
    reverse_coverages = {}
    samfile = pysam.AlignmentFile(args.input,"rb")
    logger.info(f"processing {args.input} ...")
    for contig_id in samfile.header.references:
        length = samfile.header.get_reference_length(contig_id)
        forward_coverages[contig_id] = np.zeros(length,dtype=int)
        reverse_coverages[contig_id] = np.zeros(length,dtype=int)

    strand_flip = {"forward":"reverse","reverse":"forward"}


    n = 0
    N = 0
    last_query_name = ""
    for read in samfile:
        query_name = read.query_name
        contig_id = read.reference_name     
        if query_name == last_query_name:
            if last_contig_id != contig_id:
                continue
            start2, end2 = read.reference_start,read.reference_start + read.query_alignment_length
            strand = {0:"forward",1:"reverse"}[int(read.is_reverse)]
            start = min(start1,start2)
            end = max(end1, end2)           
            N += 1
            if end - start > 1000:
                n += 1
                continue
            if N%100000 == 0:
                print(n/N)
            if ((args.strand == "forward") and (not read.is_read1)) or ((args.strand == "reverse") and (read.is_read1)) :
                strand = strand_flip[strand]                    
            if strand == "reverse":
                reverse_coverages[contig_id][start:end] += 1
            else:
                forward_coverages[contig_id][start:end] += 1  
        else:
            start1, end1 = read.reference_start,read.reference_start + read.query_alignment_length
        last_query_name = query_name
        last_contig_id = contig_id 
    scaler = None
    if args.normalize:
        scaler = 1000000/(N-n) 

    logger.info("save forward coverage ...")    
    write_bedgraph(forward_coverages, args.forward, scaler = scaler)

    logger.info("save reverse coverage ...")
    write_bedgraph(reverse_coverages, args.reverse, scaler = scaler)

    logger.info("all done .")
                
            
    

    

    

if __name__ == "__main__":
    main()
