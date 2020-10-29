#!/usr/bin/env Rscript

##load bambu
library(bambu)
args = commandArgs(trailingOnly=TRUE)

output_tag     <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore          <- strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]]
genomeseq      <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
genomeSequence <- Rsamtools::FaFile(genomeseq)
Rsamtools::indexFa(genomeseq)
annot_gtf      <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
readlist       <- args[5:length(args)]

grlist <- prepareAnnotations(annot_gtf)
se     <- bambu(reads = readlist, annotations = grlist, genomeSequence = genomeSequence, ncore = ncore, verbose = TRUE)
writeBambuOutput(se, output_tag)
