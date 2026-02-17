# Terminator-prediction-scripts
Code for "A Conserved, Protective Stem Loop Structure Irrespective of Rho Factor Dependency Enable Comparative Analysis of Bacterial Transcript 3’ Ends"

### Analysis of species with term-seq and RNA-seq data
- `run/get-coverages.sh`: Get coverage of RNA-seq reads from bam file
- `run/rho-region-finding.sh`: determine rho-dependent regions
- `run/prepare-rho-regions.sh`: define the final set of rho regions
- `run/get-folding-energy.sh`: Get folding energy of upstream sequence of transcript 3' ends, and random sequences
  
### Data augmentation of training instances
- `run/hmmsearch.sh`: Identify homologs of genes with experimentally supported terminators in at least 3 phyla
- Assign hitted genes to clusters
```{bash}
for i in {0..9};do jsub "scripts/gene-assignment.py $i";done
```
- Get candidate intervals and extract candidate sequences, and group by genes and clades for motif discovery
```{bash}
run/get-candidate-intervals.py
for i in {0..9};do "bedtools getfasta -name -s -fi otus.combined.${i}.fa -fo ${i}.fa -bed output/candidate-intervals/${i}.bed > log/get.sequence.${i}.log 2>&1";done
run/group-candidate-sequence-by-genes.py
```
- `run/run-cmfinder.sh`: run motif discovery with cmfinder
- `run/motif-filtering.sh`: discard redundant motifs, MSA with few sequences / few structures
- `run/run-cmsearch.sh`: search other gene clusters
- `run/count-hits-by-genes.sh`: count hits of each motif to other genes
- `run/get-selected-motifs.sh`: further filtering
- model calibration and search for original sequences
```{bash}
run/reformat-cm.sh
run/cailbrate-cm.sh
run/align-back.sh
scripts/get-back-align-coverage.py
```
- `run/prepare-data-1212.py`: prepare dataset for training

### Analysis of metaterm-seq data
- `snakefiles/term-seq-analysis.snakefile`: pipeline for metaterm-seq data analysis, from fastq to transcript 3' end signal in bedgraph format 
- Reads were mapped to abundant species with bowtie2
- The read 5' ends coverage (corresponds to TES signal in term-seq) was determined by a customized script based on HTSseq
- Local max positions supported by at least 10 reads were considered as experimentally supported termination sites

### Performance evaluation 
- transterm-HP
  - `scripts/run-transterm.py`: a wrapper for transterm-HP
  - `run/transterm-scan-genomes.sh`: scan genome with transterm-HP
  - `run/calibrate-transterm-refseq.sh`, `run/calibrate-transterm-GEMs.sh`: evaluate false negative rate of transterm-HP
- RDT111-OPLS-DA
  - As the model weights is not provided by the original publication (PMID:29931073), retrained one
  - `tools/mvRDT/prepare-negative-sequences.py`: simulate negative sequences by dinucleotide shuffle random genomic sequences
  - `tools/mvRDT/prepare-mvRDT-descripter.py`: a wrapper for preparing RDT111-OPLS-DA features
  - `tools/mvRDT/training.PLSDA.ensemble.py`: train the OPLS-DA model
  - `run/run-mvRDT.sh`: scanning with the OPLS-DA model
  - `/run/run-mvRDT-bg.sh`, `run/run-mvRDT-bg.sh`: scanning background sequence and evaluating false positive rate
- rhoTermPredict  
  - `run/run-rhotermpred.sh`: run rhotermpred script
  - `run/run-rhotermpred-bg.sh`,`run/get-rhotermpred-FPR-refseq.sh`: scan negative sequences and evaluate false positive rate
  - `run/get-rhotermpred-test-set-recall.sh`: performance evaluation

  
### Analysis variations of the stem loop
- `run/TES-annotation.sh`: annotate predicted 3' ends
- Training of the numeric encoder
```{bash}
#pretraining
scripts/pretrain.py --model-config config/RNA.encoder.medium.json --train-sequence dataset/pretrain/sequences.fa --train-group dataset/pretrain/clustering-table.txt --val-sequence dataset/pretrain/rRNA.fa --val-group dataset/pretrain/rRNA.txt --models models/pretrained05.weighted --performance pretrain.medium.230525.weighted.metrics.txt --device cuda:0 --weighted
# metric learning
scripts/triplet-training.py --encoder-config config/RNA.encoder.medium.json --train-sequence dataset/contrastive/rfam-full/train.fa --train-group dataset/contrastive/rfam-full/train.txt --val-sequence dataset/contrastive/rfam-full/test.fa --val-group dataset/contrastive/rfam-full/test.txt -d cuda:1 --pretrained-model models/pretrained/6.43890.pt --metrics triplet.mc.4.pc.0.1.trim.32.pos.0.8.fix.txt -m models/triplet.rfam.full.fix2 --batch-size 128 -sf 0.2 -pf 0.8 
```
- `run/term-encoding.sh`: Get the numeric embedding
  
### Analysis of putative rho binding sites
- `run/search-rho-domain.sh`: homolog search of rho proteins
- `run/extract-terminators-flanking-sequences.sh`: extract sequence flanking putative primary TES
  
### Analysis of AMR genes
- `run/run-amrfinder-GEMs.sh`: scan AMR genes with amrfinder

### Data availability
- `dataset/term-seq-sites-with-rho-annotation`: term-seq sites with rho dependency and folding energy annotation
- `dataset/term-seq-sites`: term-seq sites curated for data augmentation
- `dataset/metaterm-seq`: metaterm-seq dataset for performance evaluation
- `dataset/dRNA-seq`: dRNA-seq dataset for performance evaluation
- For a dataset too large to upload, check [this](https://cloud.tsinghua.edu.cn/d/cc891006b48b44a2b1bc/)
  - `terminators.flanked.fa.gz`: augmented terminator instances for training
  - `TES.bed.gz`: predicted stem loops in GEM representative genomes (FPR cutoff set to 0.1/KB)
  - `combined-statistics.txt`: statistics of properties in 42905 bacteria species in GEM representative genomes

  




  
  


