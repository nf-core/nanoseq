#!/usr/bin/env python

import os
import sys
import argparse

def parse_args(args=None):
    Description = 'Reformat nf-core/nanoseq samplesheet file and check its contents.'
    Epilog = """Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet file.")
    parser.add_argument('FILE_OUT', help="Output samplesheet file.")
    parser.add_argument('INPUT_PATH', help="Input path for files.")

    return parser.parse_args(args)


def print_error(error,line):
    print("ERROR: Please check samplesheet -> {}\nLine: '{}'".format(error,line.strip()))


def check_samplesheet(FileIn,FileOut,InputPath):
    HEADER = ['group','replicate','barcode','input_file','genome','transcriptome']

    ## CHECK HEADER
    fin = open(FileIn,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print("ERROR: Please check samplesheet header -> {} != {}".format(','.join(header),','.join(HEADER)))
        sys.exit(1)

    outLines, groups = [], {}
    while True:
        line = fin.readline()
        if line:
            lspl = [x.strip() for x in line.strip().split(',')]

            ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
            numCols = len([x for x in lspl[:4] if x])
            if numCols < 2:
                print_error("Please specify 'sample' entry along with either 'fastq' or 'barcode'!",line)
                sys.exit(1)

            ## ASSIGN ENTRIES TO VARIABLES
            group,replicate,barcode,input_file,genome,transcriptome = lspl

            ## CHECK GROUP ENTRIES
            if group:
                if group.find(' ') != -1:
                    print("{}: Group id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)
                else:
                    if group not in groups:
                        groups[group] = 1
                    else:
                        groups[group] += 1

            ## CHECK REPLICATE ENTRIES
            if replicate:
                if not replicate.isdigit():
                    print("{}: Replicate id not an integer!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)

            ## CHECK INPUT FILE ENTRIES
            if input_file:
                if InputPath == 'false':
                    if input_file[-9:] != '.fastq.gz' and input_file[-6:] != '.fq.gz' and input_file[-4:] != ".bam":
                        print_error("FastQ file does not have extension '.fastq.gz' or '.fq.gz' or '.bam'!",line)
                        sys.exit(1)
                    sample = input_file.split('/')[-1].split('.')[0]
                else:
                    sample = input_file
                    input_file = ''
                    if replicate:
                        sample += "REP"+replicate
            else:
                if group:
                    sample = group
                    if replicate:
                        sample += "REP"+replicate
                else:
                    print_error("Input file has not been specified!",line)
                    sys.exit(1)

            ## CHECK BARCODE ENTRIES
            if barcode:
                if not barcode.isdigit():
                    print_error("Barcode entry is not an integer!",line)
                    sys.exit(1)
                else:
                    barcode = 'barcode%s' % (barcode.zfill(2))

            ## CHECK GENOME ENTRIES
            if genome:
                if genome.find(' ') != -1:
                    print_error("Genome entry contains spaces!",line)
                    sys.exit(1)

                if len(genome.split('.')) > 1:
                    if genome[-6:] != '.fasta' and genome[-3:] != '.fa' and genome[-9:] != '.fasta.gz' and genome[-6:] != '.fa.gz':
                        print_error("Genome entry does not have extension '.fasta', '.fa', '.fasta.gz' or '.fa.gz'!",line)
                        sys.exit(1)

            ## CHECK TRANSCRIPTOME ENTRIES
            gtf = ''
            is_transcripts = '0'
            if transcriptome:

                if transcriptome.find(' ') != -1:
                    print_error("Transcriptome entry contains spaces!",line)
                    sys.exit(1)

                if transcriptome[-6:] != '.fasta' and transcriptome[-3:] != '.fa' and transcriptome[-9:] != '.fasta.gz' and transcriptome[-6:] != '.fa.gz' and transcriptome[-4:] != '.gtf' and transcriptome[-7:] != '.gtf.gz':
                    print_error("Transcriptome entry does not have extension '.fasta', '.fa', '.fasta.gz', '.fa.gz', '.gtf' or '.gtf.gz'!",line)
                    sys.exit(1)

                if transcriptome[-4:] == '.gtf' or transcriptome[-7:] == '.gtf.gz':
                    gtf = transcriptome
                    if not genome:
                        print_error("If genome isn't provided, transcriptome must be in fasta format for mapping!",line)
                        sys.exit(1)
                else:
                    is_transcripts = '1'
                    genome = transcriptome

            outLines.append([sample,group,replicate,input_file,barcode,genome,gtf,is_transcripts])
        else:
            fin.close()
            break

    ## WRITE TO FILE
    fout = open(FileOut,'w')
    fout.write(','.join(['sample','group','replicate','input_file','barcode','genome','gtf','is_transcripts']) + '\n')
    for line in outLines:
        fout.write(','.join(line) + '\n')
    fout.close()
    outfile=open("total_conditions.csv","w")
    if len(groups) > 0:
         num_samp = next(iter(groups.values()))
         if num_samp >= 3 and all(samples == num_samp for samples in groups.values()):
             outfile.write(",".join(groups)+"\n")
         else:
             outfile.write("false")
    else:
         outfile.write("false")

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN,args.FILE_OUT,args.INPUT_PATH)

if __name__ == '__main__':
    sys.exit(main())
