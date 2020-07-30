#!/usr/bin/env Rscript

##load bambu
library(bambu)

args = commandArgs(trailingOnly=TRUE)

readlist <- dir(args[1],full.names = TRUE, pattern = ".bam$")
output_tag <- args[2]
ncore <- args[3]

if (length(args) < 5) {
  stop("Please input the fullpath for the present directory and the sample sheet.", call.=FALSE)
} else {
  genomeseq <- args[5]
  genomeSequence <- Rsamtools::FaFile(genomeseq)
  Rsamtools::indexFa(genomeseq)
}
annot_gtf <- args[4]


grlist <- prepareAnnotationsFromGTF(annot_gtf)
se <- bambu(reads = readlist, annotations = grlist,genomeSequence = genomeSequence, ncore = ncore, verbose = TRUE)
writeBambuOutput(se, output_tag)

