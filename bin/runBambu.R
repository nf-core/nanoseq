#!/usr/bin/env Rscript

##install bambu if bambu is not installed
##load bambu
library(bambu)

if (!requireNamespace("BSgenome", quietly = TRUE)){
  
  install.packages("BiocManager", repos='http://cran.us.r-project.org')
  BiocManager::install("BSgenome",update = FALSE, ask= FALSE)
  
  
}


args = commandArgs(trailingOnly=TRUE)

readlist <- dir(args[1],full.names = TRUE, pattern = ".bam$")
output_tag <- args[2]
ncore <- args[3]

if (length(args) < 4) {
  stop("Please input the fullpath for the present directory and the sample sheet.", call.=FALSE)
} else if (length(args)==4) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos='http://cran.us.r-project.org')
  }
  if (!require("BSgenome.Hsapiens.NCBI.GRCh38")){
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38",update = FALSE, ask= FALSE)
    library(BSgenome.Hsapiens.NCBI.GRCh38)
  }
  genomeseq <- "BSgenome.Hsapiens.NCBI.GRCh38"  #use BSgenome if fasta file is not provided
} else {
  genomeseq <- args[5]
  genomeSequence <- Rsamtools::FaFile(genomeseq)
  Rsamtools::indexFa(genomeseq)
}
annot_gtf <- args[4]


grlist <- prepareAnnotationsFromGTF(annot_gtf)
se <- bambu(reads = readlist, annotations = grlist,genomeSequence = genomeSequence, ncore = ncore, verbose = TRUE)
writeBambuOutput(se, output_tag)

