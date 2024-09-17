#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
    ## - GENE COUNT QUANTIFICATION INPUTS FROM EITHER STRINGTIE2+FEATURECOUNTS OR BAMBU
    ## - SAMPLE INFORMATION INCLUDING CONDITIONS
    ## - THE PACKAGE BELOW NEEDS TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARY                               ##
################################################
################################################

library(DESeq2)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Please input the directory with the featureCounts results and the sample information file", call.=FALSE)
}
# default output file
outfile <- "deseq2.results.txt"

transcriptquant <- args[1]
path            <-args[2]

################################################
################################################
## FORMAT GENE COUNT QUANTIFICATION OUTPUT    ##
################################################
################################################

#create a dataframe for all samples
if (transcriptquant == "stringtie2"){
    count.matrix       <- data.frame(read.table(path, sep="\t", header=TRUE, skip = 1))
    count.matrix$Chr   <- count.matrix$Start <- count.matrix$End <- count.matrix$Length <- count.matrix$Strand <- NULL
    colnames(count.matrix)[2:length(colnames(count.matrix))] <- unlist(lapply(strsplit(colnames(count.matrix)[2:length(colnames(count.matrix))], "\\."), "[[", 1))
    count.matrix       <- aggregate(count.matrix[, -1],count.matrix["Geneid"],sum)
    countTab           <- count.matrix[, -1]
    rownames(countTab) <-count.matrix[, 1]
}
if (transcriptquant == "bambu"){
    countTab           <- data.frame(read.table(path, sep="\t", header=TRUE, row.names = 1))
    colnames(countTab) <- unlist(lapply(strsplit(colnames(countTab), "\\."), "[[", 1))
    countTab[, 1:length(colnames(countTab))] <- sapply(countTab, as.integer)
}


################################################
################################################
## READ IN SAMPLE INFORMATION (CONDITIONS)    ##
################################################
################################################

sample <- colnames(countTab)
group <- sub("(^[^-]+)_.*", "\\1", sample)
sampInfo <- data.frame(group, row.names = sample)
if (!all(rownames(sampInfo) == colnames(countTab))){
    sampInfo <- sampInfo[match(colnames(countTab), rownames(sampInfo)), ]
}

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

dds <- DESeqDataSetFromMatrix(countData = countTab, colData = sampInfo, design = ~ group)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=outfile)
