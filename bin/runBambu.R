#!/usr/bin/env Rscript

##load bambu
library(bambu)
args = commandArgs(trailingOnly=TRUE)

readlist <- strsplit(grep('--input*', args, value = TRUE), split = '=')[[1]][[2]]
readlist <- trimws(unlist(strsplit(gsub("\\[|\\]", "",readlist),",")))
output_tag <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
ncore <- strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]]

if (length(args) < 5) {
  stop("Please input the fullpath for the present directory and the sample sheet.", call.=FALSE)
} else {
  genomeseq <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
  genomeSequence <- Rsamtools::FaFile(genomeseq)
  Rsamtools::indexFa(genomeseq)
}
annot_gtf <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]


grlist <- prepareAnnotations(annot_gtf)
se <- bambu(reads = readlist, annotations = grlist,genomeSequence = genomeSequence, ncore = ncore, verbose = TRUE)
writeBambuOutput(se, output_tag)
