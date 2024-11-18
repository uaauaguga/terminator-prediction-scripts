#!/usr/bin/env python
import argparse
import numpy as np
import pysam
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('get coverage by bins')


def main():
    parser = argparse.ArgumentParser(description='get coverage by bins')
    parser.add_argument('--input','-i',type=str, required=True, help='input bam file')
    parser.add_argument('--forward','-f',type=str, required=True, help="output forward coverage file")
    parser.add_argument('--reverse','-r',type=str, required=True, help="output reverse coverage file")
    parser.add_argument('--window','-w',type=int, default=50, help="bin size to consider")
    parser.add_argument('--strand','-s',type=str, default="forward", choices=["forward","reverse"],help="libtype of the bam file")
    args = parser.parse_args()

    forward_coverages = {}
    reverse_coverages = {}
    samfile = pysam.AlignmentFile(args.input,"rb")
    logger.info(f"processing {args.input} ...")
    for contig_id in samfile.header.references:
        length = samfile.header.get_reference_length(contig_id)
        if contig_id == "AL123456.3":
            contig_id = "NC_000962.3" 
        forward_coverages[contig_id] = np.zeros(int(length/args.window)+1,dtype=int)
        reverse_coverages[contig_id] = np.zeros(int(length/args.window)+1,dtype=int)

    strand_flip = {"forward":"reverse","reverse":"forward"}
    for read in samfile:
        binidx = int((read.reference_start + read.query_alignment_length/2)/args.window)
        #binidx = int((read.reference_start + read.reference_end)/(2*args.window))
        contig_id = read.reference_name        
        if contig_id is None:          
            continue 
        if contig_id == "AL123456.3":
            contig_id = "NC_000962.3"
        strand = {0:"forward",1:"reverse"}[int(read.is_reverse)]
        if read.is_paired:        
            if ((args.strand == "forward") and (not read.is_read1)) or ((args.strand == "reverse") and (read.is_read1)) :
                strand = strand_flip[strand]                    
        else:
            if args.strand == "reverse":
                strand = strand_flip[strand]
            
        if strand == "reverse":
            reverse_coverages[contig_id][binidx] += 1
        else:
            forward_coverages[contig_id][binidx] += 1  


    ffwd = open(args.forward,"w")
    frev = open(args.reverse,"w")

    logger.info(f"saving results ...")
    for contig_id in forward_coverages:
        nbins = reverse_coverages[contig_id].shape[0]
        for i in range(nbins):
            start = i*args.window
            end = start + args.window
            print(contig_id, start, end, forward_coverages[contig_id][i],file=ffwd,sep="\t")
            print(contig_id, start, end, reverse_coverages[contig_id][i],file=frev,sep="\t")
    ffwd.close()
    frev.close()

    logger.info("all done .")

if __name__ == "__main__":
    main()
