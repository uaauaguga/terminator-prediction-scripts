#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='pick confident terminator sites')
    parser.add_argument('--input','-i',type=str, required=True, help="input reads 5' coverage in bedgraph format")
    parser.add_argument('--output','-o',type=str, required=True, help="selected putative terminator sites in bed graph format")
    parser.add_argument('--min-coverage','-mc',type=int,default=10,help="coverage required for a position to be considered as termination site")
    parser.add_argument('--min-distance','-md',type=int,default=10,help="two putative terminator site should have a distance >= this value")
    args = parser.parse_args()
    
    fout = open(args.output,"w")
    with open(args.input) as fin:
        local_max_coverage = -1
        local_max_position = -1
        last_chrom_id, last_start, last_end  = "", -1000, -1000
        for line in fin:
            fields = line.strip().split("\t")
            chrom_id, start, end, coverage = fields
            start, end, coverage = int(start), int(end), int(float(coverage))
            # ignore positions with low coverage
            if coverage < args.min_coverage:
                continue
            # update current local statistics
            if chrom_id == last_chrom_id and start - last_end < args.min_distance:
                if coverage > local_max_coverage:
                    local_max_coverage = coverage
                    local_max_position = int((start+ end)/2)
                    local_chrom_id = chrom_id
            else:
                # save current local statistics, if any
                if local_max_position >= 0:
                    print(local_chrom_id,local_max_position,local_max_position+1,".",local_max_coverage,sep="\t",file=fout)
                # start a new local statistics
                local_max_position = int((start+end)/2)
                local_max_coverage = coverage
                local_chrom_id = chrom_id
            last_chrom_id, last_start, last_end = chrom_id, start, end


if __name__ == "__main__":
    main()
