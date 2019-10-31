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
argParser.add_argument('-dm', '--demultiplex', dest="DEMULTIPLEX", help="Whether demultipexing is to be performed (default: False).",action='store_true')
argParser.add_argument('-bc', '--nobarcoding', dest="NOBARCODING", help="Whether barcode kit has been provided to Guppy (default: False).",action='store_true')
args = argParser.parse_args()

############################################
############################################
## MAIN SCRIPT
############################################
############################################

ERROR_STR = 'ERROR: Please check samplesheet'
HEADER = ['sample', 'fastq', 'barcode', 'genome']

## CHECK HEADER
fin = open(args.DESIGN_FILE_IN,'r')
header = fin.readline().strip().split(',')
if header != HEADER:
    print "{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER))
    sys.exit(1)

outLines = []
while True:
    line = fin.readline()
    if line:
        lspl = [x.strip() for x in line.strip().split(',')]
        sample,fastq,barcode,genome = lspl

        ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
        numCols = len([x for x in lspl if x])
        if numCols < 2:
            print "{}: Invalid number of columns (minimum of 2)!\nLine: '{}'".format(ERROR_STR,line.strip())
            sys.exit(1)

        if sample:
            ## CHECK SAMPLE ID HAS NO SPACES
            if sample.find(' ') != -1:
                print "{}: Sample ID contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)
        else:
            print "{}: Sample ID not specified!\nLine: '{}'".format(ERROR_STR,line.strip())
            sys.exit(1)

        if barcode:
            ## CHECK BARCODE COLUMN IS INTEGER
            if not barcode.isdigit():
                print "{}: Barcode not an integer!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)
            else:
                barcode = 'barcode%s' % (barcode.zfill(2))

        if fastq:
            ## CHECK FASTQ FILE EXTENSION
            if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                print "{}: FastQ file has incorrect extension (has to be '.fastq.gz' or 'fq.gz')!\nLine: '{}'".format(ERROR_STR,line.strip())

        if genome:
            ## CHECK GENOME HAS NO SPACES
            if genome.find(' ') != -1:
                print "{}: Genome field contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)

            ## CHECK GENOME EXTENSION
            if len(genome.split('.')) > 1:
                if genome[-6:] != '.fasta' and genome[-3:] != '.fa' and genome[-9:] != '.fasta.gz' and genome[-6:] != '.fa.gz':
                    print "{}: Genome field incorrect extension (has to be '.fasta.gz' or 'fa.gz')!\nLine: '{}'".format(ERROR_STR,line.strip())
                    sys.exit(1)

        outLines.append([sample,fastq,barcode,genome])
    else:
        fin.close()
        break

if args.NOBARCODING:
    if len(outLines) != 1:
        print "{}: Only a single-line can be specified in samplesheet without barcode information!".format(ERROR_STR)
        sys.exit(1)
    ## USE SAMPLE NAME AS BARCODE WHEN NOT DEMULTIPLEXING
    outLines[0][2] = sample

## WRITE TO FILE
fout = open(args.DESIGN_FILE_OUT,'w')
fout.write(','.join(HEADER) + '\n')
for line in outLines:
    fout.write(','.join(line) + '\n')
fout.close()
