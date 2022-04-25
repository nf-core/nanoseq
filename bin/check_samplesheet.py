#!/usr/bin/env python

import os
import sys
import errno
import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("UPDATED_PATH", help="Input path for test_nobc_nodx_rnamod")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)

class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq_1",
        second_col="fastq_2",
        single_col="single_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._single_col = single_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        assert len(row[self._first_col]) > 0, "At least the first FASTQ file is required."
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            assert (
                Path(row[self._first_col]).suffixes == Path(row[self._second_col]).suffixes
            ), "FASTQ pairs must have the same file extensions."
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FORMATS), (
            f"The FASTQ file has an unrecognized extension: {filename}\n"
            f"It should be one of: {', '.join(self.VALID_FORMATS)}"
        )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename the sample if more than one sample,
        FASTQ file combination exists.

        """
        assert len(self._seen) == len(self.modified), "The pair of sample name and FASTQ must be unique."
        if len({pair[0] for pair in self._seen}) < len(self._seen):
            counts = Counter(pair[0] for pair in self._seen)
            seen = Counter()
            for row in self.modified:
                sample = row[self._sample_col]
                seen[sample] += 1
                if counts[sample] > 1:
                    row[self._sample_col] = f"{sample}_T{seen[sample]}"


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

def check_samplesheet(file_in, updated_path, file_out):
    """
    peek = handle.read(2048)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical(f"The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    handle.seek(0)
    return dialect

    group,replicate,barcode,input_file,genome,transcriptome
    MCF7,1,,MCF7_directcDNA_replicate1.fastq.gz,genome.fa,
    MCF7,2,,MCF7_directcDNA_replicate3.fastq.gz,genome.fa,genome.gtf
    K562,1,,K562_directcDNA_replicate1.fastq.gz,genome.fa,
    K562,2,,K562_directcDNA_replicate4.fastq.gz,,transcripts.fa
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

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
            nanopolish_fast5 = ''
            if input_file:
                if input_file.find(" ") != -1:
                    print_error("Input file contains spaces!", 'Line', line)
                if input_file.endswith(".fastq.gz"):
                    input_extensions.append("*.fastq.gz")
                elif input_file.endswith(".fq.gz"):
                    input_extensions.append("*.fq.gz")
                elif input_file.endswith(".bam"):
                    input_extensions.append("*.bam")
                else:
                    if updated_path != "not_changed":
                        input_file='/'.join([updated_path,input_file.split("/")[-1]])
                    list_dir         = os.listdir(input_file)
                    nanopolish_fast5 = input_file
                    if not (all(fname.endswith('.fast5') for fname in list_dir)):
                        if "fast5" in list_dir and "fastq" in list_dir:
                            nanopolish_fast5 = input_file+'/fast5'
                            ## CHECK FAST5 DIRECTORY
                            if not (all(fname.endswith('.fast5') for fname in os.listdir(nanopolish_fast5))):
                                print_error('fast5 directory contains non-fast5 files.')
                            ## CHECK PROVIDED BASECALLED FASTQ
                            fastq_path       = input_file+'/fastq'
                            basecalled_fastq = [fn for fn in os.listdir(fastq_path) if fn.endswith(".fastq.gz") or fn.endswith(".fq.gz") ]
                            print(basecalled_fastq)
                            if len(basecalled_fastq) != 1:
                                print_error('Please input one basecalled fastq per sample.')
                            else:
                                input_file   = fastq_path+'/'+basecalled_fastq[0]
                                if not basecalled_fastq[0].endswith(".fastq.gz"):
                                    if not basecalled_fastq[0].endswith(".fq.gz"):
                                        print_error('basecalled fastq input does not end with ".fastq.gz" or ".fq.gz"')
                        else:
                            print_error('path does not end with ".fastq.gz", ".fq.gz", or ".bam" and is not an existing directory with correct fast5 and/or fastq inputs.')

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

            ## Create sample mapping dictionary = {group: {replicate : [ barcode, input_file, genome, gtf, is_transcripts, nanopolish_fast5 ]}}
            sample_info = [barcode, input_file, genome, gtf, is_transcripts, nanopolish_fast5]
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

            fout.write(",".join(['sample', 'barcode', 'input_file', 'genome', 'gtf', 'is_transcripts', 'nanopolish_fast5']) + "\n")
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
    check_samplesheet(args.FILE_IN, args.UPDATED_PATH, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
