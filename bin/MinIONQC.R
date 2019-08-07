#!/usr/bin/env Rscript

# MinIONQC version 1.3.5
# Copyright (C) 2017 onwards Robert Lanfear
#
# For license see https://github.com/roblanf/minion_qc/blob/master/LICENSE


# supress warnings
options(warn=-1)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))




# option parsing #

parser <- OptionParser()

parser <- add_option(parser, 
                     opt_str = c("-i", "--input"), 
                     type = "character",
                     dest = 'input.file',
                     help="Input file or directory (required). Either a full path to a sequence_summary.txt file, or a full path to a directory containing one or more such files. In the latter case the directory is searched recursively."
                     )

parser <- add_option(parser, 
                     opt_str = c("-o", "--outputdirectory"), 
                     type = "character",
                     dest = 'output.dir',
                     default=NA,
                     help="Output directory (optional, default is the same as the input directory). If a single sequencing_summary.txt file is passed as input, then the output directory will contain just the plots associated with that file. If a directory containing more than one sequencing_summary.txt files is passed as input, then the plots will be put into sub-directories that have the same names as the parent directories of each sequencing_summary.txt file"
                     )

parser <- add_option(parser, 
                     opt_str = c("-q", "--qscore_cutoff"), 
                     type="double", 
                     default=7.0,
                     dest = 'q',
                     help="The cutoff value for the mean Q score of a read (default 7). Used to create separate plots for reads above and below this threshold"
                     )

parser <- add_option(parser, 
                     opt_str = c("-p", "--processors"), 
                     type="integer", 
                     default=1,
                     dest = 'cores',
                     help="Number of processors to use for the anlaysis (default 1). Only helps when you are analysing more than one sequencing_summary.txt file at a time"
                     )

parser <- add_option(parser, 
                     opt_str = c("-s", "--smallfigures"), 
                     type = "logical",
                     default = FALSE,
                     dest = 'smallfig',
                     help="TRUE or FALSE (the default). When true, MinIONQC will output smaller figures, e.g. suitable for publications or presentations. The default is to produce larger figures optimised for display on screen. Some figures just require small text, and cannot be effectively resized."
                     )

parser <- add_option(parser, 
                     opt_str = c("-c", "--combined-only"), 
                     type = "logical",
                     default = FALSE,
                     dest = 'combined_only',
                     help="TRUE or FALSE (the default). When true, MinIONQC will only produce the combined report, it will not produce individual reports for each flowcell."
                     )

opt = parse_args(parser)

if(exists("test.file") == FALSE){
    test.file = c(1, 2, 3) # dummy variable
}

if (length(opt$input.file)==1) {
    input.file = opt$input.file
} else if (length(test.file)==1) {
    input.file = test.file # specifically for testing the script
    flog.info(paste("Using test file", test.file))
} else {
    stop("Input file parameter must be supplied via -i or --input. See script usage (--help) or readme for help: https://github.com/roblanf/minion_qc")
}

q = opt$q
cores = opt$cores
smallfig = opt$smallfig
combined_only = opt$combined_only

p1m = 1.0 


# this is how we label the reads at least as good as q
q_title = paste("Q>=", q, sep="")


# build the map for R9.5
p1 = data.frame(channel=33:64, row=rep(1:4, each=8), col=rep(1:8, 4))
p2 = data.frame(channel=481:512, row=rep(5:8, each=8), col=rep(1:8, 4))
p3 = data.frame(channel=417:448, row=rep(9:12, each=8), col=rep(1:8, 4))
p4 = data.frame(channel=353:384, row=rep(13:16, each=8), col=rep(1:8, 4))
p5 = data.frame(channel=289:320, row=rep(17:20, each=8), col=rep(1:8, 4))
p6 = data.frame(channel=225:256, row=rep(21:24, each=8), col=rep(1:8, 4))
p7 = data.frame(channel=161:192, row=rep(25:28, each=8), col=rep(1:8, 4))
p8 = data.frame(channel=97:128, row=rep(29:32, each=8), col=rep(1:8, 4))

q1 = data.frame(channel=1:32, row=rep(1:4, each=8), col=rep(16:9, 4))
q2 = data.frame(channel=449:480, row=rep(5:8, each=8), col=rep(16:9, 4))
q3 = data.frame(channel=385:416, row=rep(9:12, each=8), col=rep(16:9, 4))
q4 = data.frame(channel=321:352, row=rep(13:16, each=8), col=rep(16:9, 4))
q5 = data.frame(channel=257:288, row=rep(17:20, each=8), col=rep(16:9, 4))
q6 = data.frame(channel=193:224, row=rep(21:24, each=8), col=rep(16:9, 4))
q7 = data.frame(channel=129:160, row=rep(25:28, each=8), col=rep(16:9, 4))
q8 = data.frame(channel=65:96, row=rep(29:32, each=8), col=rep(16:9, 4))

map = rbind(p1, p2, p3, p4, p5, p6, p7, p8, q1, q2, q3, q4, q5, q6, q7, q8)

add_cols <- function(d, min.q){
    # take a sequencing sumamry file (d), and a minimum Q value you are interested in (min.q)
    # return the same data frame with the following columns added 
        # cumulative.bases and cumulative.bases.time
        # hour of run
        # reads.per.hour
    
    d = subset(d, mean_qscore_template >= min.q)
    
    if(nrow(d)==0){ 
        flog.error(paste("There are no reads with a mean Q score higher than your cutoff of ", min.q, ". Please choose a lower cutoff and try again.", sep = ""))
        quit()
    }
    
    if(max(d$channel)<=512){
        d = merge(d, map, by="channel")
    }else{
        # thanks to Matt Loose. Code adapted from: https://github.com/mattloose/flowcellvis/blob/master/flowcellgif.py
        block = floor((d$channel-1)/250)
        remainder = (d$channel-1)%%250
        d$row = floor(remainder/10) + 1 # +1 because R is not zero indexed
        d$col = remainder%%10 + block*10 + 1 # +1 because R is not zero indexed 
    }
    
    d = d[with(d, order(as.numeric(start_time))), ] # sort by start time
    d$cumulative.bases.time = cumsum(as.numeric(d$sequence_length_template))

    d = d[with(d, order(-sequence_length_template)), ] # sort by read length
    d$cumulative.bases = cumsum(as.numeric(d$sequence_length_template))
    d$hour = d$start_time %/% 3600
    
    # add the reads generated for each hour
    reads.per.hour = as.data.frame(table(d$hour))
    names(reads.per.hour) = c("hour", "reads_per_hour")
    reads.per.hour$hour = as.numeric(as.character(reads.per.hour$hour))
    d = merge(d, reads.per.hour, by = c("hour"))    
    return(d)
}


load_summary <- function(filepath, min.q){
    # load a sequencing summary and add some info
    # min.q is a vector of length 2 defining 2 levels of min.q to have
    # by default the lowest value is -Inf, i.e. includes all reads. The 
    # other value in min.q is set by the user at the command line
    d = read_tsv(filepath, col_types = cols_only(channel = 'i', 
                                                num_events_template = 'i', 
                                                sequence_length_template = 'i', 
                                                mean_qscore_template = 'n',
                                                sequence_length_2d = 'i',
                                                mean_qscore_2d = 'n',
                                                start_time = 'n',
                                                calibration_strand_genome_template = 'c'))
    
    if(max(d$channel)<=512){
        flog.info("MinION flowcell detected")
    }else{
        flog.info("PromethION flowcell detected")
    }
    
    # remove the control sequence from directRNA runs
    if("calibration_strand_genome_template" %in% names(d)){
        d = subset(d, calibration_strand_genome_template != "YHR174W")
    }

    if("sequence_length_2d" %in% names(d)){
        # it's a 1D2 or 2D run
        d$sequence_length_template = as.numeric(as.character(d$sequence_length_2d))
        d$mean_qscore_template = as.numeric(as.character(d$mean_qscore_2d))
        d$num_events_template = NA
        d$start_time = as.numeric(as.character(d$start_time))
        
    }else{
        d$sequence_length_template = as.numeric(as.character(d$sequence_length_template))
        d$mean_qscore_template = as.numeric(as.character(d$mean_qscore_template))
        d$num_events_template = as.numeric(as.character(d$num_events_template))
        d$start_time = as.numeric(as.character(d$start_time))
    }
        
    d$events_per_base = d$num_events_template/d$sequence_length_template

    flowcell = basename(dirname(filepath))
    
    # add columns for all the reads
    d1 = add_cols(d, min.q[1])
    d1$Q_cutoff = "All reads"
    
    # add columns for just the reads that pass the user Q threshold
    d2 = add_cols(d, min.q[2])
    d2$Q_cutoff = q_title
    
    # bind those two together into one data frame
    d = as.data.frame(rbindlist(list(d1, d2)))

    # name the flowcell (useful for analyses with >1 flowcell)
    d$flowcell = flowcell
    
    # make sure this is a factor
    d$Q_cutoff = as.factor(d$Q_cutoff)
    
    keep = c("hour", "start_time", "channel", "sequence_length_template", "mean_qscore_template", "row", "col", "cumulative.bases", "cumulative.bases.time", "reads_per_hour", "Q_cutoff", "flowcell", "events_per_base")
    dk = d[, which(names(d) %in% keep)]
        
    return(dk)
}

reads.gt <- function(d, len){
    # return the number of reads in data frame d
    # that are at least as long as length len
    return(length(which(d$sequence_length_template>=len)))
}

bases.gt <- function(d, len){
    # return the number of bases contained in reads from 
    # data frame d
    # that are at least as long as length len
    reads = subset(d, sequence_length_template >= len)
    return(sum(as.numeric(reads$sequence_length_template)))
}

log10_minor_break = function (...){
    # function to add minor breaks to a log10 graph
    # hat-tip: https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
    function(x) {
        minx         = floor(min(log10(x), na.rm=T))-1;
        maxx         = ceiling(max(log10(x), na.rm=T))+1;
        n_major      = maxx-minx+1;
        major_breaks = seq(minx, maxx, by=1)
        minor_breaks = 
            rep(log10(seq(1, 9, by=1)), times = n_major)+
            rep(major_breaks, each = 9)
        return(10^(minor_breaks))
    }
}


log10_major_break = function (...){
    # function to add major breaks to a log10 graph
    # hat-tip: https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
    function(x) {
        minx         = floor(min(log10(x), na.rm=T))-1;
        maxx         = ceiling(max(log10(x), na.rm=T))+1;
        n_major      = maxx-minx+1;
        major_breaks = seq(minx, maxx, by=1)
        return(10^(major_breaks))
    }
}



binSearch <- function(min, max, df, t = 100000) {
    # binary search algorithm, thanks to https://stackoverflow.com/questions/46292438/optimising-a-calculation-on-every-cumulative-subset-of-a-vector-in-r/46303384#46303384
    # the aim is to return the number of reads in a dataset (df)
    # that comprise the largest subset of reads with an N50 of t
    # we use this to calculte the number of 'ultra long' reads
    # which are defined as those with N50 > 100KB
    mid = floor(mean(c(min, max)))
    if (mid == min) {
        if (df$sequence_length_template[min(which(df$cumulative.bases>df$cumulative.bases[min]/2))] < t) {
            return(min - 1)
        } else {
            return(max - 1)
        }
    }
    
    n = df$sequence_length_template[min(which(df$cumulative.bases>df$cumulative.bases[mid]/2))]
    if (n >= t) {
        return(binSearch(mid, max, df))
    } else {
        return(binSearch(min, mid, df))
    }
}


summary.stats <- function(d, Q_cutoff="All reads"){
    # Calculate summary stats for a single value of min.q
    
    rows = which(as.character(d$Q_cutoff)==Q_cutoff)
    d = d[rows,]
    d = d[with(d, order(-sequence_length_template)), ] # sort by read length, just in case
    
    total.bases = sum(as.numeric(d$sequence_length_template))
    total.reads = nrow(d)
    N50.length = d$sequence_length_template[min(which(d$cumulative.bases > (total.bases/2)))]
    mean.length = round(mean(as.numeric(d$sequence_length_template)), digits = 1)
    median.length = round(median(as.numeric(d$sequence_length_template)), digits = 1)
    max.length = max(as.numeric(d$sequence_length_template))
    mean.q = round(mean(d$mean_qscore_template), digits = 1)
    median.q = round(median(d$mean_qscore_template), digits = 1)
    
    #calculate ultra-long reads and bases (max amount of data with N50>100KB)
    ultra.reads = binSearch(1, nrow(d), d, t = 100000)    
    if(ultra.reads>=1){
        ultra.gigabases = sum(as.numeric(d$sequence_length_template[1:ultra.reads]))/1000000000
    }else{
        ultra.gigabases = 0
    }
        
    reads = list(
                reads.gt(d, 10000), 
                reads.gt(d, 20000), 
                reads.gt(d, 50000),
                reads.gt(d, 100000),
                reads.gt(d, 200000),
                reads.gt(d, 500000),
                reads.gt(d, 1000000),
                ultra.reads)
    names(reads) = c(">10kb", ">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m", "ultralong")

    bases = list(
                bases.gt(d, 10000)/1000000000, 
                bases.gt(d, 20000)/1000000000, 
                bases.gt(d, 50000)/1000000000,
                bases.gt(d, 100000)/1000000000,
                bases.gt(d, 200000)/1000000000,
                bases.gt(d, 500000)/1000000000,
                bases.gt(d, 1000000)/1000000000,
                ultra.gigabases)
    names(bases) = c(">10kb", ">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m", "ultralong")
    
    return(list('total.gigabases' = total.bases/1000000000,
                'total.reads' = total.reads,
                'N50.length' = N50.length, 
                'mean.length' = mean.length, 
                'median.length' = median.length,
                'max.length' = max.length,
                'mean.q' = mean.q,
                'median.q' = median.q,
                'reads' = reads,
                'gigabases' = bases
                ))
}

channel.summary <- function(d){
    # calculate summaries of what happened in each of the channels 
    # of a flowcell
    
    a = ddply(d, .(channel), 
              summarize, 
              total.bases = sum(sequence_length_template), 
              total.reads = sum(which(sequence_length_template>=0)), 
              mean.read.length = mean(sequence_length_template), 
              median.read.length = median(sequence_length_template),
              row = mean(row),
              col = mean(col))
    b = melt(a, id.vars = c("channel", "row", "col"))
    return(b)    
}


single.flowcell <- function(input.file, output.dir, q=7, base.dir = NA){
    # wrapper function to analyse data from a single flowcell
    # input.file is a sequencing_summary.txt file from a 1D run
    # output.dir is the output directory into which to write results
    # q is the cutoff used for Q values, set by the user
    # base.dir is the base directory if and only if the user supplied a base directory
    # we use base.dir to name flowcells in a sensible way
    flog.info(paste("Loading input file:", input.file))
    d = load_summary(input.file, min.q=c(-Inf, q))

    flowcell = unique(d$flowcell)

    
    # output goes with the sequencing summary file unless otherwise specified
    if(is.na(opt$output.dir)){ 
        output.dir = file.path(dirname(input.file))
    } else {
        # the user supplied an output dir
        output.dir = file.path(opt$output.dir, flowcell)
    }
    
    flog.info(paste(sep = "", flowcell, ": creating output directory:", output.dir))
    dir.create(output.dir, recursive = TRUE)
    out.txt = file.path(output.dir, "summary.yaml")
    
    flog.info(paste(sep = "", flowcell, ": summarising input file for flowcell"))
    all.reads.summary = summary.stats(d, Q_cutoff = "All reads")
    q10.reads.summary = summary.stats(d, Q_cutoff = q_title)
    
    summary = list("input file" = input.file,
                   "All reads" = all.reads.summary,
                   cutoff = q10.reads.summary,
                   "notes" = 'ultralong reads refers to the largest set of reads with N50>100KB')
    
    names(summary)[3] = q_title
    
    write(as.yaml(summary), out.txt)
    
    muxes = seq(from = 0, to = max(d$hour), by = 8)

    # set up variable sizes
    if(smallfig == TRUE){ p1m = 0.5 }else{ p1m = 1.0 }
    if(smallfig == TRUE){ p2m = 0.6 }else{ p2m = 1.0 }

    
    # make plots
    flog.info(paste(sep = "", flowcell, ": plotting length histogram"))
    p1 = ggplot(d, aes(x = sequence_length_template, fill = Q_cutoff)) + 
        geom_histogram(bins = 300) + 
        scale_x_log10(minor_breaks=log10_minor_break(), breaks = log10_major_break()) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Read length (bases)") +
        ylab("Number of reads") +
        guides(fill=FALSE) + scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75)
    suppressMessages(ggsave(filename = file.path(output.dir, "length_histogram.png"), width = p1m*960/75, height = p1m*960/75, plot = p1)) #

    flog.info(paste(sep = "", flowcell, ": plotting mean Q score histogram"))
    p2 = ggplot(d, aes(x = mean_qscore_template, fill = Q_cutoff)) + 
        geom_histogram(bins = 300) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Mean Q score of read") +
        ylab("Number of reads") +
        guides(fill=FALSE) + scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75)
    suppressMessages(ggsave(filename = file.path(output.dir, "q_histogram.png"), width = p1m*960/75, height = p1m*960/75, plot = p2)) #
    
    if(max(d$channel)<=512){
        # only do this for minion, not promethion
	    flog.info(paste(sep = "", flowcell, ": plotting flowcell overview"))
        p3 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x=start_time/3600, y=sequence_length_template, colour = mean_qscore_template)) + 
            geom_point(size=1.5, alpha=0.35) + 
            scale_colour_viridis() + 
            labs(colour='Q')  + 
            scale_y_log10() + 
            facet_grid(row~col) +
            theme(panel.spacing = unit(0.5, "lines")) +
            xlab("Hours into run") +
            ylab("Read length") +
            theme(text = element_text(size = 40), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), legend.text=element_text(size=18), legend.title=element_text(size=24))
    
            suppressMessages(ggsave(filename = file.path(output.dir, "flowcell_overview.png"), width = 2000/75, height = 1920/75, plot = p3))
    }

    flog.info(paste(sep = "", flowcell, ": plotting flowcell yield over time"))
    p5 = ggplot(d, aes(x=start_time/3600, y=cumulative.bases.time/1000000000, colour = Q_cutoff)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        geom_line(size = 1) + 
        xlab("Hours into run") +
        ylab("Total yield in gigabases") +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        theme(text = element_text(size = 15))
    suppressMessages(ggsave(filename = file.path(output.dir, "yield_over_time.png"), width = p1m*960/75, height = p1m*480/75, plot = p5)) #


    flog.info(paste(sep = "", flowcell, ": plotting flowcell yield by read length"))
    p6 = ggplot(d, aes(x=sequence_length_template, y=cumulative.bases/1000000000, colour = Q_cutoff)) + 
        geom_line(size = 1) + 
        xlab("Minimum read length (bases)") +
        ylab("Total yield in gigabases") +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        theme(text = element_text(size = 15))
    xmax = max(d$sequence_length_template[which(d$cumulative.bases > 0.01 * max(d$cumulative.bases))])
    p6 = p6 + scale_x_continuous(limits = c(0, xmax))
    suppressMessages(ggsave(filename = file.path(output.dir, "yield_by_length.png"), width = p1m*960/75, height = p1m*480/75, plot = p6)) #
    
    flog.info(paste(sep = "", flowcell, ": plotting sequence length over time"))
    p7 = ggplot(d, aes(x=start_time/3600, y=sequence_length_template, colour = Q_cutoff, group = Q_cutoff)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean read length (bases)") + 
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        ylim(0, NA)
    suppressMessages(ggsave(filename = file.path(output.dir, "length_by_hour.png"), width = p1m*960/75, height = p1m*480/75, plot = p7)) #
    
    flog.info(paste(sep = "", flowcell, ": plotting Q score over time"))
    p8 = ggplot(d, aes(x=start_time/3600, y=mean_qscore_template, colour = Q_cutoff, group = Q_cutoff)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean Q score") + 
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        ylim(0, NA)
    suppressMessages(ggsave(filename = file.path(output.dir, "q_by_hour.png"), width = p1m*960/75, height = p1m*480/75, plot = p8)) #
    
    flog.info(paste(sep = "", flowcell, ": plotting reads per hour"))
    f = d[c("hour", "reads_per_hour", "Q_cutoff")]
    f = f[!duplicated(f),]
    g = subset(f, Q_cutoff=="All reads")
    h = subset(f, Q_cutoff==q_title)
    max = max(f$hour)
    # all of this is just to fill in hours with no reads recorded
    all = 0:max
    add.g = all[which(all %in% g$hour == FALSE)]
    if(length(add.g)>0){
        add.g = data.frame(hour = add.g, reads_per_hour = 0, Q_cutoff = "All reads")
        g = rbind(g, add.g)
    }
    add.h = all[which(all %in% h$hour == FALSE)]
    if(length(add.h)>0){
        add.h = data.frame(hour = add.h, reads_per_hour = 0, Q_cutoff = q_title)
        h = rbind(h, add.h)
    }
    i = rbind(g, h)
    i$Q_cutoff = as.character(i$Q_cutoff)
    i$Q_cutoff[which(i$Q_cutoff==q_title)] = paste("Q>=", q, sep="")
    p9 = ggplot(i, aes(x=hour, y=reads_per_hour, colour = Q_cutoff, group = Q_cutoff)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_point() +
        geom_line() +
        xlab("Hours into run") + 
        ylab("Number of reads per hour") + 
        ylim(0, NA) + 
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads"))
    suppressMessages(ggsave(filename = file.path(output.dir, "reads_per_hour.png"), width = p1m*960/75, height = p1m*480/75, plot = p9)) #
    
    if(max(d$channel)<=512){
        # minion
        flog.info(paste(sep = "", flowcell, ": plotting read length vs. q score scatterplot"))
        p10 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + 
            geom_point(alpha=0.05, size = 0.4) + 
            scale_x_log10(minor_breaks=log10_minor_break(), breaks = log10_major_break()) + 
            labs(colour='Events per base\n(log scale)\n')  + 
            theme(text = element_text(size = 15)) +
            xlab("Read length (bases)") +
            ylab("Mean Q score of read")
    }else{
        # promethion
        p10 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + 
            geom_bin2d() + 
            scale_x_log10(minor_breaks=log10_minor_break(), breaks = log10_major_break()) + 
            theme(text = element_text(size = 15)) +
            scale_fill_viridis() +
            xlab("Read length (bases)") +
            ylab("Mean Q score of read")
    }
    
    
    if(max(d$events_per_base, na.rm=T)>0){
        # a catch for 1D2 runs which don't have events per base
        p10 = p10 + scale_colour_viridis(trans = "log", labels = scientific, option = 'inferno') 
    }
    # we keep it a bit wider, because the legend takes up a fair bit of the plot space
    suppressMessages(ggsave(filename = file.path(output.dir, "length_vs_q.png"), width = p2m*960/75, height = p1m*960/75, plot = p10)) #
    
    flog.info(paste(sep = "", flowcell, ": plotting flowcell channels summary histograms"))
    c = channel.summary(subset(d, Q_cutoff=="All reads"))
    c10 = channel.summary(subset(d, Q_cutoff==q_title))
    c$Q_cutoff = "All reads"
    c10$Q_cutoff = q_title
    cc = rbind(c, c10)
    cc$variable = as.character(cc$variable)
    cc$variable[which(cc$variable=="total.bases")] = "Number of bases per channel"
    cc$variable[which(cc$variable=="total.reads")] = "Number of reads per channel"
    cc$variable[which(cc$variable=="mean.read.length")] = "Mean read length per channel"
    cc$variable[which(cc$variable=="median.read.length")] = "Median read length per channel"
    p11 = ggplot(cc, aes(x = value, fill = Q_cutoff)) + geom_histogram(bins = 30) + 
        facet_grid(Q_cutoff~variable, scales = "free_x") + 
        theme(text = element_text(size = 15), axis.text.x = element_text(angle = 60, hjust = 1)) +
        guides(fill=FALSE) +
        scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75) +
        guides(fill=FALSE)
    suppressMessages(ggsave(filename = file.path(output.dir, "channel_summary.png"), width = 960/75, height = 480/75, plot = p11)) 
    

    flog.info(paste(sep = "", flowcell, ": plotting physical overview of output per channel"))
    if(max(d$channel)<=512){
        # minion
        cols = 2
    }else{
        # promethion
        cols = 1
    }


    p12 = ggplot(subset(cc, variable == "Number of bases per channel"), aes(x = as.numeric(col), y = as.numeric(row))) + 
        geom_tile(aes(fill = value/1000000000), colour="white", size=0.25) +
        facet_wrap(~Q_cutoff, ncol = cols) + 
        theme(text = element_text(size = 15), 
                plot.background=element_blank(), 
                panel.border=element_blank(),
                panel.background = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()) +
        scale_fill_viridis(name = "GB/channel") +
        scale_y_continuous(trans = "reverse", expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0)) +
        coord_fixed() +
        labs(x="channel column",y="channel row")

    if(max(d$channel)<=512){
        # minion
        suppressMessages(ggsave(filename = file.path(output.dir, "gb_per_channel_overview.png"), width = 960/150, height = 480/75, plot = p12)) 
    }else{
        # promethion
        suppressMessages(ggsave(filename = file.path(output.dir, "gb_per_channel_overview.png"), width = 960/75, height = 480/75, plot = p12))         
    }
    
    return(d)
}

combined.flowcell <- function(d, output.dir, q=8){
    # function to analyse combined data from multiple flowcells
    # useful for getting an overall impression of the combined data

    flog.info("Creating output directory")
    out.txt = file.path(output.dir, "summary.yaml")
    
    # write summaries
    flog.info(paste("Summarising combined data from all flowcells, saving to:", out.txt))

    # tidy up and remove added stuff
    drops = c("cumulative.bases", "hour", "reads.per.hour")
    d = d[ , !(names(d) %in% drops)]
    
    d1 = subset(d, Q_cutoff == "All reads")    
    d1 = d1[with(d1, order(-sequence_length_template)), ] # sort by read length
    d1$cumulative.bases = cumsum(as.numeric(d1$sequence_length_template))
        
    d2 = subset(d, Q_cutoff == q_title)    
    d2 = d2[with(d2, order(-sequence_length_template)), ] # sort by read length
    d2$cumulative.bases = cumsum(as.numeric(d2$sequence_length_template))

    d1$Q_cutoff = as.factor(d1$Q_cutoff)
    d2$Q_cutoff = as.factor(d2$Q_cutoff)
    
    all.reads.summary = summary.stats(d1, Q_cutoff = "All reads")
    q10.reads.summary = summary.stats(d2, Q_cutoff = q_title)
    
    summary = list("input file" = input.file,
                   "All reads" = all.reads.summary,
                   cutoff = q10.reads.summary,
                   "notes" = 'ultralong reads refers to the largest set of reads with N50>100KB')
    
    names(summary)[3] = q_title
    
    write(as.yaml(summary), out.txt)

    d = rbind(d1, d2)
    d$Q_cutoff = as.factor(d$Q_cutoff)
    d1 = 0
    d2 = 0
        
    # set up variable sizes
    if(smallfig == TRUE){ p1m = 0.5 }else{ p1m = 1.0 }

    # make plots
    flog.info("Plotting combined length histogram")
    p1 = ggplot(d, aes(x = sequence_length_template, fill = Q_cutoff)) + 
        geom_histogram(bins = 300) + 
        scale_x_log10(minor_breaks=log10_minor_break(), breaks = log10_major_break()) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Read length (bases)") +
        ylab("Number of reads") +
        guides(fill=FALSE) + scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75)
    suppressMessages(ggsave(filename = file.path(output.dir, "combined_length_histogram.png"), width = p1m*960/75, height = p1m*960/75, plot = p1))
    
    flog.info("Plotting combined mean Q score histogram")
    p2 = ggplot(d, aes(x = mean_qscore_template, fill = Q_cutoff)) + 
        geom_histogram(bins = 300) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Mean Q score of read") +
        ylab("Number of reads") +
        guides(fill=FALSE) + scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75)
    suppressMessages(ggsave(filename = file.path(output.dir, "combined_q_histogram.png"), width = p1m*960/75, height = p1m*960/75, plot = p2))
        
    flog.info("Plotting combined yield by length")
    p4 = ggplot(d, aes(x=sequence_length_template, y=cumulative.bases/1000000000, colour = Q_cutoff)) + 
        geom_line(size = 1) + 
        xlab("Minimum read length (bases)") +
        ylab("Total yield in gigabases") +
        scale_colour_viridis(discrete = TRUE, begin = 0.25, end = 0.75, guide = guide_legend(title = "Reads")) +
        theme(text = element_text(size = 15))
    xmax = max(d$sequence_length_template[which(d$cumulative.bases > 0.01 * max(d$cumulative.bases))])
    p4 = p4 + scale_x_continuous(limits = c(0, xmax))
    
    suppressMessages(ggsave(filename = file.path(output.dir, "combined_yield_by_length.png"), width = p1m*960/75, height = p1m*480/75, plot = p4))
    
}


multi.flowcell = function(input.file, output.base, q){
    # wrapper function to allow parallelisation of single-flowcell 
    # analyses when >1 flowcell is analysed in one run
    
    d = single.flowcell(input.file, output.dir, q)
    return(d)
    
}


multi.plots = function(dm, output.dir){
    # function to plot data from multiple flowcells,
    # where the data is not combined (as in combined.flowcell() )
    # but instead just uses multiple lines on each plot.
    
    muxes = seq(from = 0, to = max(dm$hour), by = 8)
    
    # set up variable sizes
    if(smallfig == TRUE){ p1m = 0.5 }else{ p1m = 1.0 }

    # make plots
    flog.info("Plotting length distributions")
    p1 = ggplot(dm, aes(x = sequence_length_template)) + 
        geom_line(stat="density", aes(colour = flowcell), size = 1) +
        scale_x_log10(minor_breaks=log10_minor_break(), breaks = log10_major_break()) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Read length (bases)") +
        ylab("Density")
    suppressMessages(ggsave(filename = file.path(output.dir, "length_distributions.png"), width = p1m*960/75, height = p1m*960/75, plot = p1)) #
    
    flog.info("Plotting mean Q score distributions")
    p2 = ggplot(dm, aes(x = mean_qscore_template)) + 
        geom_line(stat="density", aes(colour = flowcell), size = 1) +
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Mean Q score of read") +
        ylab("Density")
    suppressMessages(ggsave(filename = file.path(output.dir, "q_distributions.png"), width = p1m*960/75, height = p1m*960/75, plot = p2)) #


    flog.info("Plotting flowcell yield over time")
    p5 = ggplot(dm, aes(x=start_time/3600, y=cumulative.bases.time/1000000000, colour = flowcell)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        geom_line(size = 1) + 
        xlab("Hours into run") +
        ylab("Total yield in gigabases") +
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") +
        theme(text = element_text(size = 15))
    suppressMessages(ggsave(filename = file.path(output.dir, "yield_over_time.png"), width = p1m*960/75, height = p1m*960/75, plot = p5)) #

        
    flog.info("Plotting flowcell yield by length")
    p6 = ggplot(dm, aes(x=sequence_length_template, y=cumulative.bases/1000000000, colour = flowcell)) + 
        geom_line(size = 1) + 
        xlab("Minimum read length (bases)") +
        ylab("Total yield in gigabases") +
        theme(text = element_text(size = 15)) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y")
    xmax = max(dm$sequence_length_template[which(dm$cumulative.bases > 0.01 * max(dm$cumulative.bases))])
    p6 = p6 + scale_x_continuous(limits = c(0, xmax))
    suppressMessages(ggsave(filename = file.path(output.dir, "yield_by_length.png"), width = p1m*960/75, height = p1m*960/75, plot = p6)) #
    

    flog.info("Plotting sequence length over time")
    p7 = ggplot(dm, aes(x=start_time/3600, y=sequence_length_template, colour = flowcell)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean read length (bases)") + 
        ylim(0, NA) +
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y")
    suppressMessages(ggsave(filename = file.path(output.dir, "length_by_hour.png"), width = p1m*960/75, height = p1m*960/75, plot = p7))
    

    flog.info("Plotting Q score over time")
    p8 = ggplot(dm, aes(x=start_time/3600, y=mean_qscore_template, colour = flowcell)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean Q score") + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y")
    suppressMessages(ggsave(filename = file.path(output.dir, "q_by_hour.png"), width = p1m*960/75, height = p1m*960/75, plot = p8))    
    
}


# Choose how to act depending on whether we have a single input file or mulitple input files

if(file_test("-f", input.file)==TRUE & length(test.file)>1){
    # if it's an existing file (not a folder) just run one analysis
    d = single.flowcell(input.file, output.dir, q)
}else if(file_test("-d", input.file)==TRUE & length(test.file)>1){
    # it's a directory, recursively analyse all sequencing_summary.txt files
  
    # get a list of all sequencing_summary.txt files, recursively
    summaries = list.files(path = input.file, pattern = "sequencing_summary.txt", recursive = TRUE, full.names = TRUE)

    flog.info("")
    flog.info("**** Analysing the following files ****")
    flog.info(summaries)

    
    # if the user passes a directory with only one sequencing_summary.txt file...
    if(length(summaries) == 1){
        d = single.flowcell(summaries[1], output.dir, q)
        flog.info('**** Analysis complete ****')
    }else{
        # analyse each one and keep the returns in a list
        if(combined_only == FALSE){
            results = mclapply(summaries, multi.flowcell, output.dir, q, mc.cores = cores)
        }else{
            results = mclapply(summaries, load_summary, min.q = c(-Inf, q), mc.cores = cores)
        }
        
        # rbind that list
        flog.info('**** Analysing data from all flowcells combined ****')
        dm = as.data.frame(rbindlist(results))

        # now do the single plot on ALL the output
        if(is.na(opt$output.dir)){ 
            combined.output = file.path(input.file, "combinedQC")
        } else {
            # the user supplied an output dir
            combined.output = file.path(opt$output.dir, "combinedQC")
        }
        
        flog.info(paste("Plots from the combined output will be saved in", combined.output))
        dir.create(combined.output, recursive = TRUE)
        combined.flowcell(dm, combined.output, q)
        multi.plots(dm, combined.output)
        flog.info('**** Analysis complete ****')
        flog.info('If you use MinIONQC in your published work, please cite:')
        flog.info('R Lanfear, M Schalamun, D Kainer, W Wang, B Schwessinger (2018). MinIONQC: fast and simple quality control for MinION sequencing data, Bioinformatics, bty654')
        flog.info('https://doi.org/10.1093/bioinformatics/bty654')
    }
    
}else{
    #WTF
    flog.warn(paste("Couldn't find a sequencing summary file in your input which was: ", 
        input.file, 
        "\nThe input must be either a sequencing_summary.txt file, or a directory containing one or more such files"))
}

