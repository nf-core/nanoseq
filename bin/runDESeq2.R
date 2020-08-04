#!/usr/bin/env Rscript

library(DESeq2)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Please input the directory with the featureCounts results and the sample information file", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  args[4] = "DESeq2out.txt"
}
#DeSeq2
transcriptquant <- args[1]
path<-args[2]
#create a dataframe for all samples
if (transcriptquant == "stringtie2"){
  #count_files<- grep(list.files(path), pattern='tx_', inv=TRUE, value=TRUE)
  count.matrix <- data.frame(read.table(path,sep="\t",header=TRUE, skip = 1))
  count.matrix$Chr <- count.matrix$Start <- count.matrix$End <- count.matrix$Length <- count.matrix$Strand <- NULL
  colnames(count.matrix)[2:length(colnames(count.matrix))] <- unlist(lapply(strsplit(colnames(count.matrix)[2:length(colnames(count.matrix))],"\\."),"[[",1))
  count.matrix <- aggregate(count.matrix[,-1],count.matrix["Geneid"],sum)
  countTab <- count.matrix[,-1]
  rownames(countTab)<-count.matrix[,1]
}
if (transcriptquant == "bambu"){
   countTab <- data.frame(read.table(path,sep="\t",header=TRUE))
   colnames(countTab) <- unlist(lapply(strsplit(colnames(countTab),"\\."),"[[",1))
}

#sampInfo <- read.csv("~/Downloads/nanorna-bam-master/two_conditions.csv",row.names = 1)
sampInfo<-read.csv(args[3],row.names=1)
if (!all(rownames(sampInfo) == colnames(countTab))){
  sampInfo <- sampInfo[match(colnames(countTab), rownames(sampInfo)),]
}
dds <- DESeqDataSetFromMatrix(countData = countTab,colData = sampInfo,design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
#register(MulticoreParam(6))
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=args[4])
