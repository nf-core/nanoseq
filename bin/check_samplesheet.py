#!/usr/bin/env python

# TODO nf-core: Update the script to check the samplesheet
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/nanoseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    group,replicate,barcode,input_file,genome,transcriptome
    MCF7,1,,MCF7_directcDNA_replicate1.fastq.gz,genome.fa,
    MCF7,2,,MCF7_directcDNA_replicate3.fastq.gz,genome.fa,genome.gtf
    K562,1,,K562_directcDNA_replicate1.fastq.gz,genome.fa,
    K562,2,,K562_directcDNA_replicate4.fastq.gz,,transcripts.fa
    """

    input_extensions = []
    sample_info_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 3
        HEADER = ['group', 'replicate', 'barcode', 'input_file', 'genome', 'transcriptome']
        header = fin.readline().strip().split(",")
        if header[:len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip() for x in line.strip().split(",")]

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error("Invalid number of columns (minimum = {})!".format(len(HEADER)), 'Line', line)

            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error("Invalid number of populated columns (minimum = {})!".format(MIN_COLS), 'Line', line)

            ## Check group name entries
            group, replicate, barcode, input_file, genome, transcriptome = lspl[:len(HEADER)]
            if group:
                if group.find(" ") != -1:
                    print_error("Group entry contains spaces!", 'Line', line)
            else:
                print_error("Group entry has not been specified!", 'Line', line)

            ## Check replicate entry is integer
            if replicate:
                if not replicate.isdigit():
                    print_error("Replicate id not an integer!", 'Line', line)
            else:
                print_error("Replicate id not specified!", 'Line', line)
            replicate = int(replicate)

            ## Check barcode entry
            if barcode:
                if not barcode.isdigit():
                    print_error("Barcode entry is not an integer!", 'Line', line)
                else:
                    barcode = 'barcode%s' % (barcode.zfill(2))

            ## Check input file extension
            if input_file:
                if input_file.find(" ") != -1:
                    print_error("Input file contains spaces!", 'Line', line)
                if not input_file.endswith(".fastq.gz") and not input_file.endswith(".fq.gz") and not input_file.endswith(".bam"):
                    print_error("Input file does not have extension '.fastq.gz', '.fq.gz' or '.bam'!", 'Line', line)
                if input_file.endswith(".fastq.gz"):
                    input_extensions.append("*.fastq.gz")
                elif input_file.endswith(".fq.gz"):
                    input_extensions.append("*.fq.gz")
                elif input_file.endswith(".bam"):
                    input_extensions.append("*.bam")

            ## Check genome entries
            if genome:
                if genome.find(' ') != -1:
                    print_error("Genome entry contains spaces!",'Line', line)
                if len(genome.split('.')) > 1:
                    if genome[-6:] != '.fasta' and genome[-3:] != '.fa' and genome[-9:] != '.fasta.gz' and genome[-6:] != '.fa.gz':
                        print_error("Genome entry does not have extension '.fasta', '.fa', '.fasta.gz' or '.fa.gz'!",'Line', line)

            ## Check transcriptome entries
            gtf = ''
            is_transcripts = '0'
            if transcriptome:
                if transcriptome.find(' ') != -1:
                    print_error("Transcriptome entry contains spaces!",'Line',line)
                if transcriptome[-6:] != '.fasta' and transcriptome[-3:] != '.fa' and transcriptome[-9:] != '.fasta.gz' and transcriptome[-6:] != '.fa.gz' and transcriptome[-4:] != '.gtf' and transcriptome[-7:] != '.gtf.gz':
                    print_error("Transcriptome entry does not have extension '.fasta', '.fa', '.fasta.gz', '.fa.gz', '.gtf' or '.gtf.gz'!", 'Line', line)
                if transcriptome[-4:] == '.gtf' or transcriptome[-7:] == '.gtf.gz':
                    gtf = transcriptome
                    if not genome:
                        print_error("If genome isn't provided, transcriptome must be in fasta format for mapping!", 'Line', line)
                else:
                    is_transcripts = '1'
                    genome = transcriptome

            ## Create sample mapping dictionary = {group: {replicate : [ barcode, input_file, genome, gtf, is_transcripts ]}}
            sample_info = [barcode, input_file, genome, gtf, is_transcripts]
            if group not in sample_info_dict:
                sample_info_dict[group] = {}
            if replicate not in sample_info_dict[group]:
                sample_info_dict[group][replicate] = sample_info
            else:
                print_error("Same replicate id provided multiple times!", 'Line', line)

    ## Check all input files have the same extension
    if len(set(input_extensions)) > 1:
        print_error("All input files must have the same extension!", 'Multiple extensions found', ', '.join(set(input_extensions)))

    ## Write validated samplesheet with appropriate columns
    if len(sample_info_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:

            fout.write(",".join(['sample', 'barcode', 'input_file', 'genome', 'gtf', 'is_transcripts']) + "\n")
            for sample in sorted(sample_info_dict.keys()):

                ## Check that replicate ids are in format 1..<NUM_REPS>
                uniq_rep_ids = set(sample_info_dict[sample].keys())
                if len(uniq_rep_ids) != max(uniq_rep_ids):
                    print_error("Replicate ids must start with 1..<num_replicates>!", 'Group', sample)

                ### Write to file
                for replicate in sorted(sample_info_dict[sample].keys()):
                    sample_id = "{}_R{}".format(sample,replicate)
                    fout.write(','.join([sample_id] + sample_info_dict[sample][replicate]) + '\n')


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
