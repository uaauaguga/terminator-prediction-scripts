#!/usr/bin/env Rscript
library(ape)
library(phylolm)

tree <- read.tree("cyanobacteria-with-complete-genome-and-traits.marker120.asm-id.nwk")
metadata <- read.table("metadata-complete-cleaned.txt",sep="\t",header=TRUE,row.names=1)
metadata$X01.MARINE <- 1-metadata$X01.NON.MARINE
res <- phyloglm(X01.MARINE ~ RDT.fraction, metadata, tree, method = "logistic_IG10",start.alpha=0.1)
summary(res)
res <- phyloglm(X02.NO.N.FIXATION ~ RDT.fraction, metadata, tree, method = "logistic_IG10",start.alpha=0.1)
summary(res)
res <- phyloglm(X03.UNICELLULAR ~ RDT.fraction, metadata, tree, method = "logistic_IG10",start.alpha=0.1)
summary(res)
