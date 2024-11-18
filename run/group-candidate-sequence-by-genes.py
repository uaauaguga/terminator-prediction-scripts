#!/usr/bin/env python
import pandas as  pd

def main():
    taxonomy = pd.read_csv("/apps/home/lulab_jinyunfan/qhsky1/metagenome-ncRNA/otu-metadata/otu_taxonomy.tsv",sep="\t",index_col=0)
    otu_ids = taxonomy[taxonomy["t_domain"] == "Bacteria"].index
    otu_ids = set(list(otu_ids))
    groupped_sequences = {}
    for i in range(10):
        path = f"output/candidate-sequences/chunks/{i}.fa"
        print(f"Load {i} ...")
        with open(path) as f:
            for line in f:
                # >OTU-45000:NZ_JQMO01000001.1_5:0012::OTU-45000:NZ_JQMO01000001.1:6363-6507(-)
                seq_id = line.strip()[1:]
                gene_id = seq_id.split("::")[0]
                fields = gene_id.split(":")
                otu_id = fields[0]
                gene_id = fields[-1]
                sequence = next(f).strip()
                if otu_id not in otu_ids:
                    continue
                if len(sequence) <= 75:
                    continue
                if gene_id not in groupped_sequences:
                    groupped_sequences[gene_id] = []
                groupped_sequences[gene_id].append((seq_id,sequence))

    for gene_id in groupped_sequences:
        if len(groupped_sequences[gene_id]) < 1000:
            continue
        path = f"output/candidate-sequences/by-genes/{gene_id}.fa" 
        print(f"Save cluster {gene_id} ...")
        with open(path,"w") as f:
            for seq_id, sequence in groupped_sequences[gene_id]:
                f.write(f">{seq_id}\n")
                f.write(f"{sequence}\n")


if __name__ == "__main__":
    main()
