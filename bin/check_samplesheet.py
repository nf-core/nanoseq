#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on August 28th 2019 to check nf-core/nanoseq design file
#######################################################################
#######################################################################

import os
import sys
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Reformat nf-core/nanoseq design file and check its contents.'
Epilog = """Example usage: python check_samplesheet.py <DESIGN_FILE_IN> <DESIGN_FILE_OUT>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('DESIGN_FILE_IN', help="Input design file.")
argParser.add_argument('DESIGN_FILE_OUT', help="Output design file.")

## OPTIONAL PARAMETERS
argParser.add_argument('-sd', '--skip_demultiplexing', dest="SKIP_DEMULTIPLEXING", help="Whether demultipexing is to be performed (default: False).",action='store_true')
args = argParser.parse_args()

############################################
############################################
## MAIN SCRIPT
############################################
############################################

ERROR_STR = 'ERROR: Please check samplesheet'
HEADER = ['sample', 'fastq', 'barcode', 'genome', 'transcriptome']

## CHECK HEADER
fin = open(args.DESIGN_FILE_IN,'r')
header = fin.readline().strip().split(',')
if header != HEADER:
    print("{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER)))
    sys.exit(1)

outLines = []
while True:
    line = fin.readline()
    if line:
        lspl = [x.strip() for x in line.strip().split(',')]
        sample,fastq,barcode,genome,transcriptome = lspl

        ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
        numCols = len([x for x in lspl if x])
        if numCols < 2:
            print("{}: Invalid number of columns (minimum of 2)!\nLine: '{}'".format(ERROR_STR,line.strip()))
            sys.exit(1)

        ## CHECK SAMPLE ID ENTRIES
        if sample:
            if sample.find(' ') != -1:
                print("{}: Sample ID contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
        else:
            print("{}: Sample ID not specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
            sys.exit(1)

        ## CHECK BARCODE ENTRIES
        if barcode:
            if not barcode.isdigit():
                print("{}: Barcode not an integer!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            else:
                barcode = 'barcode%s' % (barcode.zfill(2))

        ## CHECK FASTQ ENTRIES
        if fastq:
            if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                print("{}: FastQ file has incorrect extension (has to be '.fastq.gz' or '.fq.gz')!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

        ## CHECK GENOME ENTRIES
        if genome:
            if genome.find(' ') != -1:
                print("{}: Genome field contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

            if len(genome.split('.')) > 1:
                if genome[-6:] != '.fasta' and genome[-3:] != '.fa' and genome[-9:] != '.fasta.gz' and genome[-6:] != '.fa.gz':
                    print("{}: Genome field incorrect extension (has to be '.fasta', '.fa', '.fasta.gz' or '.fa.gz')!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)

        ## CHECK TRANSCRIPTOME ENTRIES
        gtf = ''
        is_transcripts = '0'
        if transcriptome:

            if transcriptome.find(' ') != -1:
                print("{}: Transcriptome field contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

            if transcriptome[-6:] != '.fasta' and transcriptome[-3:] != '.fa' and transcriptome[-9:] != '.fasta.gz' and transcriptome[-6:] != '.fa.gz' and transcriptome[-4:] != '.gtf':
                print("{}: Transcriptome field incorrect extension (has to be '.fasta', '.fa', '.fasta.gz', '.fa.gz' or '.gtf')!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

            if transcriptome[-4:] == '.gtf':
                gtf = transcriptome
                if not genome:
                    print("{}: If genome isnt provided, transcriptome must be in fasta format for mapping!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)
            else:
                is_transcripts = '1'
                genome = transcriptome
                
        outLines.append([sample,fastq,barcode,genome,gtf,is_transcripts])
    else:
        fin.close()
        break

if args.SKIP_DEMULTIPLEXING:
    if len(outLines) != 1:
        print("{}: Only a single-line can be specified in samplesheet without barcode information!".format(ERROR_STR))
        sys.exit(1)
    ## USE SAMPLE NAME AS BARCODE WHEN NOT DEMULTIPLEXING
    outLines[0][2] = outLines[0][0]

## WRITE TO FILE
fout = open(args.DESIGN_FILE_OUT,'w')
fout.write(','.join(['sample', 'fastq', 'barcode', 'genome', 'gtf', 'is_transcripts']) + '\n')
for line in outLines:
    fout.write(','.join(line) + '\n')
fout.close()
