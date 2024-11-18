#!/bin/bash
scripts/rhoterm-region-finder.py --treatments output/B.B/coverage/SRR22959760.+.txt,output/B.B/coverage/SRR22959761.+.txt --controls output/B.B/coverage/SRR22959762.+.txt,output/B.B/coverage/SRR22959763.+.txt --strand + --output output/B.B/rho-regions.+.txt
scripts/rhoterm-region-finder.py --treatments output/B.B/coverage/SRR22959760.-.txt,output/B.B/coverage/SRR22959761.-.txt --controls output/B.B/coverage/SRR22959762.-.txt,output/B.B/coverage/SRR22959763.-.txt --strand - --output output/B.B/rho-regions.-.txt


scripts/rhoterm-region-finder.py --treatments output/M.T/coverage/rnaseq_rhoDUC_3h_r1.+.txt,output/M.T/coverage/rnaseq_rhoDUC_3h_r2.+.txt --controls output/M.T/coverage/rnaseq_rhoDUC_0h_r1.+.txt,output/M.T/coverage/rnaseq_rhoDUC_0h_r2.+.txt --output output/M.T/rho-regions.3.vs.0.+.txt  --strand + 
scripts/rhoterm-region-finder.py --treatments output/M.T/coverage/rnaseq_rhoDUC_3h_r1.-.txt,output/M.T/coverage/rnaseq_rhoDUC_3h_r2.-.txt --controls output/M.T/coverage/rnaseq_rhoDUC_0h_r1.-.txt,output/M.T/coverage/rnaseq_rhoDUC_0h_r2.-.txt --output output/M.T/rho-regions.3.vs.0.-.txt  --strand - 

# finally use this
scripts/rhoterm-region-finder.py --treatments output/M.T/coverage/rnaseq_rhoDUC_6h_r1.+.txt,output/M.T/coverage/rnaseq_rhoDUC_6h_r2.+.txt --controls output/M.T/coverage/rnaseq_rhoDUC_0h_r1.+.txt,output/M.T/coverage/rnaseq_rhoDUC_0h_r2.+.txt --output output/M.T/rho-regions.6.vs.0.+.txt  --strand + 
scripts/rhoterm-region-finder.py --treatments output/M.T/coverage/rnaseq_rhoDUC_6h_r1.-.txt,output/M.T/coverage/rnaseq_rhoDUC_6h_r2.-.txt --controls output/M.T/coverage/rnaseq_rhoDUC_0h_r1.-.txt,output/M.T/coverage/rnaseq_rhoDUC_0h_r2.-.txt --output output/M.T/rho-regions.6.vs.0.-.txt  --strand - 

scripts/rhoterm-region-finder.py --treatments output/M.T/coverage/rnaseq_rhoDUC_45h_r1.+.txt,output/M.T/coverage/rnaseq_rhoDUC_45h_r2.+.txt --controls output/M.T/coverage/rnaseq_rhoDUC_0h_r1.+.txt,output/M.T/coverage/rnaseq_rhoDUC_0h_r2.+.txt --output output/M.T/rho-regions.45.vs.0.+.txt  --strand + 
scripts/rhoterm-region-finder.py --treatments output/M.T/coverage/rnaseq_rhoDUC_45h_r1.-.txt,output/M.T/coverage/rnaseq_rhoDUC_45h_r2.-.txt --controls output/M.T/coverage/rnaseq_rhoDUC_0h_r1.-.txt,output/M.T/coverage/rnaseq_rhoDUC_0h_r2.-.txt --output output/M.T/rho-regions.45.vs.0.-.txt  --strand -

scripts/rhoterm-region-finder.py --treatments output/B.S/coverage/SRR17789022.+.txt,output/B.S/coverage/SRR17789023.+.txt --controls output/B.S/coverage/SRR17789026.+.txt,output/B.S/coverage/SRR17789027.+.txt  --output output/B.S/rho-regions.+.txt  --strand + 
scripts/rhoterm-region-finder.py --treatments output/B.S/coverage/SRR17789022.-.txt,output/B.S/coverage/SRR17789023.-.txt --controls  output/B.S/coverage/SRR17789026.-.txt,output/B.S/coverage/SRR17789027.-.txt  --output output/B.S/rho-regions.-.txt  --strand - 


scripts/rhoterm-region-finder.py --treatments output/E.C/coverage/SRR609954.+.txt,output/E.C/coverage/SRR609955.+.txt  --controls output/E.C/coverage/SRR609952.+.txt,output/E.C/coverage/SRR609953.+.txt --output output/E.C/rho-regions.+.txt  --strand +
scripts/rhoterm-region-finder.py --treatments output/E.C/coverage/SRR609954.-.txt,output/E.C/coverage/SRR609955.-.txt  --controls output/E.C/coverage/SRR609952.-.txt,output/E.C/coverage/SRR609953.-.txt --output output/E.C/rho-regions.-.txt  --strand -
