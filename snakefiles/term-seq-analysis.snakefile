sample_ids = open("sample_ids.txt").read().strip().split("\n") 
rule all:
    input:
        abundance = expand("output/metaphlan/bam/{sample_id}.bam",sample_id=sample_ids),
        bigwig = expand('output/coverage/{sample_id}.{strand}.bigwig',sample_id=sample_ids,strand=["+","-"]),
        bigwig5 = expand('output/5p-coverage/{sample_id}.{strand}.bigwig',sample_id=sample_ids,strand=["+","-"]),
        bam = expand("output/bam/{sample_id}.bam",sample_id=sample_ids),
        sites = expand('output/sites/{sample_id}.{strand}.bed',sample_id=sample_ids,strand=["+","-"])
rule metaphlan:
    input:
        fastq = 'data/fastq/{sample_id}.fastq.gz',
    output:
        abundance = "output/metaphlan/abundance/{sample_id}.txt",
        bt2output = "output/metaphlan/hits/{sample_id}.txt.gz",
        bam = "output/metaphlan/bam/{sample_id}.bam"
    threads: 4
    params:
        bt2db = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases"
    log: 'output/log/metaphlan/{sample_id}.txt'
    shell:
        """
        metaphlan {input.fastq} --min_ab 0.001  --input_type fastq  --bowtie2db {params.bt2db} --bowtie2out output/metaphlan/hits/{wildcards.sample_id}.txt --output_file {output.abundance} --nproc {threads} --tmp_dir tmp --samout output/metaphlan/bam/{wildcards.sample_id}.sam
        gzip output/metaphlan/hits/{wildcards.sample_id}.txt
        samtools view -Sb output/metaphlan/bam/{wildcards.sample_id}.sam  > {output.bam}
        rm output/metaphlan/bam/{wildcards.sample_id}.sam
        """


rule mapping:
    input:
        fastq = 'data/fastq/{sample_id}.fastq.gz'
    output:
        bam = "output/bam/{sample_id}.bam",
        unmapped = "output/unmapped/{sample_id}.fastq.gz"
    threads: 6
    log: "log/mapping/{sample_id}.log"
    shell:
        """
        bowtie2 -p {threads} -U {input.fastq} -x genome/bt2-index/genomes --sensitive-local --no-unal --un-gz {output.unmapped} 2> {log} | samtools sort --output-fmt BAM --threads {threads} -o {output.bam} 
        samtools index {output.bam}
        """

rule get_coverage:
    input:
        bam = "output/bam/{sample_id}.bam",
        chromsize = "genome/fasta/chrom.size"
    output:
        fwd_bedgraph = 'output/coverage/{sample_id}.+.bedgraph',
        fwd_bigwig   = 'output/coverage/{sample_id}.+.bigwig',
        rev_bedgraph = 'output/coverage/{sample_id}.-.bedgraph',
        rev_bigwig   = 'output/coverage/{sample_id}.-.bigwig'
    shell:
        """
        bedtools genomecov -du -strand + -ibam {input.bam} -bg -split | LC_ALL=C sort -k1,1 -k2,2n > {output.fwd_bedgraph}
        bedGraphToBigWig {output.fwd_bedgraph} {input.chromsize} {output.fwd_bigwig} 
        bedtools genomecov -du -strand - -ibam {input.bam} -bg -split | LC_ALL=C sort -k1,1 -k2,2n > {output.rev_bedgraph}
        bedGraphToBigWig {output.rev_bedgraph} {input.chromsize} {output.rev_bigwig}
        """

rule get_5p_coverage:
    input:
        bam = "output/bam/{sample_id}.bam",
        chromsize = "genome/fasta/chrom.size"
    output:
        fwd_bigwig = 'output/5p-coverage/{sample_id}.+.bigwig',
        rev_bigwig = 'output/5p-coverage/{sample_id}.-.bigwig'
        fw_bg = 'output/5p-coverage/{sample_id}.sorted.+.bedgraph',
        rev_bg = 'output/5p-coverage/{sample_id}.sorted.-.bedgraph',
    log: "log/five-prime-coverage/{sample_id}.log"
    shell:
        """
        scripts/get-5p-coverage.py -i {input.bam} -f output/5p-coverage/{wildcards.sample_id}.+.bedgraph -r output/5p-coverage/{wildcards.sample_id}.-.bedgraph > {log} 2>&1
        cat output/5p-coverage/{wildcards.sample_id}.+.bedgraph  | grep -v 'type=bedGraph' | LC_ALL=C sort -k1,1 -k2,2n > output/5p-coverage/{wildcards.sample_id}.sorted.+.bedgraph
        cat output/5p-coverage/{wildcards.sample_id}.-.bedgraph  | grep -v 'type=bedGraph' | LC_ALL=C sort -k1,1 -k2,2n > output/5p-coverage/{wildcards.sample_id}.sorted.-.bedgraph
        bedGraphToBigWig output/5p-coverage/{wildcards.sample_id}.sorted.+.bedgraph {input.chromsize} {output.fwd_bigwig}
        bedGraphToBigWig output/5p-coverage/{wildcards.sample_id}.sorted.-.bedgraph {input.chromsize} {output.rev_bigwig}
        """

rule pick_tes:
    input:
        fwd_bg = 'output/5p-coverage/{sample_id}.sorted.+.bedgraph',
        rev_bg = 'output/5p-coverage/{sample_id}.sorted.-.bedgraph'
    output:
        fwd_tes = 'output/sites/{sample_id}.+.bed',
        rev_tes = 'output/sites/{sample_id}.-.bed'
    shell:
        """
        scripts/pick-term-sites.py -i {input.fwd_bg} -o {output.fwd_tes}
        scripts/pick-term-sites.py -i {input.rev_bg} -o {output.rev_tes}
        """
