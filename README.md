# Terminator-prediction-scripts
Code for "A Conserved, Protective Stem Loop Structure Irrespective of Rho Factor Dependency Enable Comparative Analysis of Bacterial Transcript 3’ Ends"

- Analysis of species with term-seq and RNA-seq data
  - `run/get-coverages.sh`: Get coverage of RNA-seq reads from bam file
  - `run/rho-region-finding.sh`: determine rho dependent regions
  - `run/prepare-rho-regions.sh`: define the final set of rho regions
  - `run/get-folding-energy.sh`: Get folding energy of upstream sequence of transcript 3' ends, and random sequences
  
- Data augmentation of traning instances
  
  
- Performance evaluation  

  
- Analysis variations of the stem loop
  - `run/TES-annotation.sh`: annotate predicted 3' ends
  - Traning of the numeric encoder
```{bash}
#pretraining
scripts/pretrain.py --model-config config/RNA.encoder.medium.json --train-sequence dataset/pretrain/sequences.fa --train-group dataset/pretrain/clustering-table.txt --val-sequence dataset/pretrain/rRNA.fa --val-group dataset/pretrain/rRNA.txt --models models/pretrained05.weighted --performance pretrain.medium.230525.weighted.metrics.txt --device cuda:0 --weighted
# metric learning
scripts/triplet-training.py --encoder-config config/RNA.encoder.medium.json --train-sequence dataset/contrastive/rfam-full/train.fa --train-group dataset/contrastive/rfam-full/train.txt --val-sequence dataset/contrastive/rfam-full/test.fa --val-group dataset/contrastive/rfam-full/test.txt -d cuda:1 --pretrained-model models/pretrained/6.43890.pt --metrics triplet.mc.4.pc.0.1.trim.32.pos.0.8.fix.txt -m models/triplet.rfam.full.fix2 --batch-size 128 -sf 0.2 -pf 0.8 
```
  - `run/term-encoding.sh`: Get the numeric embedding
  
- Analysis putative rho binding sites
  - `run/extract-terminators-flanking-sequences.sh`: extract sequence flanking putative primary TES
  
  
- Analysis of AMR genes

  
  


